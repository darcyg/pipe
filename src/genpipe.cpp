#include "genpipe.h"

GenPipe::GenPipe(Options* opt, Pipeline* p){
    mOpt = opt;
    mPipe = p;
}

GenPipe::~GenPipe(){
}

void GenPipe::genPrelibPipe(){
    std::ifstream fr(mOpt->clOpt.sample_list);
    std::string sampleNo, flowCell, libName, read1, read2, barcodeConf, tmpStr;
    std::getline(fr, tmpStr);
    std::istringstream iss;
    int count = 0;
    while(std::getline(fr, tmpStr)){
        iss.clear();
        iss.str(tmpStr);
        iss >> sampleNo >> flowCell >> libName >> read1 >> read2 >> barcodeConf;
        // generate jobs
        GenJob* genJob = new GenJob(mOpt);
        // generate spliter job
        Job* jSpliter = new Job("spliter", mOpt->ioOpt.spl_dir, libName, 1);
        genJob->setLib(read1, read2);
        genJob->genSpliterJob(barcodeConf, jSpliter);
        // push to tasks
        Task* task = new Task(1);
        task->addJob(jSpliter, 0);
        std::string sjmStatus(mOpt->ioOpt.sjm_dir + "/" + libName + "_prelib.sjm.status");
        std::map<std::string, std::string> jMap;
        Job::getStatus(jMap, sjmStatus);
        for(auto& e: task->joblist){
            for(auto& f: e){
                if(std::find(mOpt->clOpt.ana_marker.cbegin(), mOpt->clOpt.ana_marker.cend(), f->stage_marker) == mOpt->clOpt.ana_marker.cend()){
                    f->status.second = "done";
                }
                if(mOpt->clOpt.update && jMap.find(f->name.second) != jMap.end()){
                    f->status.second = jMap[f->name.second];
                }
                if(!mOpt->clOpt.queue.empty()){
                    f->queue.second = mOpt->clOpt.queue;
                }
            }
        }
        std::string sjmFile(mOpt->ioOpt.sjm_dir + "/" + libName + "_prelib.sjm");
        std::string sjmLog(mOpt->ioOpt.sjm_dir + "/" + libName + "_prelib.log");
        RunTask* runTask = new RunTask();
        runTask->sjmCMD = mOpt->ioOpt.bin_dir + "/sjm -i -l " + sjmLog + " " + sjmFile;
        runTask->goodMarkFile = mOpt->ioOpt.log_dir + "/" + libName + ".Prelib.SUCCESS";
        runTask->failMarkFile = mOpt->ioOpt.log_dir + "/" + libName + ".Prelib.FAIL";
        runTask->logFile = sjmLog;
        mPipe->addRunFile(runTask, count, 0);
        ++count;
        std::ofstream fw(sjmFile);
        fw << task;
        fw.close();
    }
    fr.close();
}

void GenPipe::genAnalibTask(){
    std::ifstream fr1(mOpt->clOpt.sample_list);
    std::string sampleNo, flowCell, libName, read1, read2, barcodeConf;
    std::string subLib, subRead1, subRead2, tmpStr;
    std::getline(fr1, tmpStr);
    std::istringstream iss1, iss2;
    int count = 0;
    while(std::getline(fr1, tmpStr)){
        iss1.clear();
        iss1.str(tmpStr);
        iss1 >> sampleNo >> flowCell >> libName >> read1 >> read2 >> barcodeConf;
        std::ifstream fr2(barcodeConf);
        while(std::getline(fr2, tmpStr)){
            iss2.clear();
            iss2.str(tmpStr);
            iss2 >> subLib;
            subRead1 = mOpt->ioOpt.spl_dir + "/" + subLib + ".R1.fq";
            subRead2 = mOpt->ioOpt.spl_dir + "/" + subLib + ".R2.fq";
            GenJob* genJob = new GenJob(mOpt);
            Task* task = new Task(8);
            // fqtool
            Job* jFqtool = new Job("fqtool", mOpt->ioOpt.cut_dir, subLib, 2);
            genJob->setLib(subRead1, subRead2);
            genJob->genFqtoolJob(jFqtool);
            task->addJob(jFqtool, 0);
            // filter
            Job* jFiltdb = new Job("filter", mOpt->ioOpt.fil_dir, subLib, 3);
            genJob->setLib(jFqtool->o1, jFqtool->o2);
            genJob->genFilterJob(jFiltdb);
            task->addJob(jFiltdb, 1); 
            // seqtk
            Job* jSeqtk = new Job("seqtk", mOpt->ioOpt.dfq_dir, subLib, 4);
            jSeqtk->optPre = sampleNo;
            genJob->setLib(jFiltdb->o1, jFiltdb->o2);
            genJob->genSeqtkJob(jSeqtk);
            task->addJob(jSeqtk, 2);
            // bwa mem
            Job* jAln = new Job("bwa", mOpt->ioOpt.aln_dir, subLib, 5);
            genJob->setLib(jSeqtk->o1, jSeqtk->o2);
            genJob->genAlignJob(jAln);
            task->addJob(jAln, 3);
            // mkdup
            Job* jMkdup = new Job("duplexer", mOpt->ioOpt.mkd_dir, subLib, 6);
            genJob->setBam(jAln->o1);
            genJob->genMkdupJob(jMkdup);
            task->addJob(jMkdup, 4);
            // bamqc
            Job* jBamqc = new Job("bamqc", mOpt->ioOpt.bqc_dir, subLib, 7);
            genJob->setBam(jMkdup->o1);
            genJob->genBamqcJob(jBamqc);
            task->addJob(jBamqc, 5);
            // fusion
            Job* jFusion = new Job("fusionmap", mOpt->ioOpt.fus_dir, subLib, 8);
            genJob->setLib(jSeqtk->o1, jSeqtk->o2);
            genJob->genFusionJob(jFusion);
            task->addJob(jFusion, 3);
            // express
            Job* jExpress = new Job("kallisto", mOpt->ioOpt.exp_dir, subLib, 9);
            genJob->setLib(jSeqtk->o1, jSeqtk->o2);
            genJob->genExpressJob(jExpress);
            task->addJob(jExpress, 3);
            // report
            Job* jReport = new Job("genrpt", mOpt->ioOpt.rep_dir, subLib, 10);
            genJob->genReportJob(jReport);
            task->addJob(jReport, 6);
            // cleanup
            if(!mOpt->clOpt.noclean){
                Job* jClean = new Job("cleanup", mOpt->ioOpt.log_dir, subLib, 11);
                genJob->genCleanupJob(jClean);
                task->addJob(jClean, 7);
            }
            // updating status
            std::string sjmStatus = mOpt->ioOpt.sjm_dir + "/" + subLib + "_analib.sjm.status";
            std::map<std::string, std::string> jMap;
            Job::getStatus(jMap, sjmStatus);
            for(auto& e: task->joblist){
                for(auto& f: e){
                    if(std::find(mOpt->clOpt.ana_marker.cbegin(), mOpt->clOpt.ana_marker.cend(), f->stage_marker) == mOpt->clOpt.ana_marker.cend()){
                        f->status.second = "done";
                    }
                    if(mOpt->clOpt.update && jMap.find(f->name.second) != jMap.end()){
                        f->status.second = jMap[f->name.second];
                    }
                    if(!mOpt->clOpt.queue.empty()){
                        f->queue.second = mOpt->clOpt.queue;
                    }
                }
            }
            std::string sjmFile(mOpt->ioOpt.sjm_dir + "/" + subLib + "_analib.sjm");
            std::string sjmLog(mOpt->ioOpt.sjm_dir + "/" + subLib + "_analib.log");
            RunTask* runTask = new RunTask();
            runTask->sjmCMD = mOpt->ioOpt.bin_dir + "/sjm -i -l " + sjmLog + " " + sjmFile;
            runTask->goodMarkFile = mOpt->ioOpt.log_dir + "/" + subLib + ".Analysis.SUCCESS";
            runTask->failMarkFile = mOpt->ioOpt.log_dir + "/" + subLib + ".Analysis.FAIL";
            runTask->logFile = sjmLog;
            mPipe->addRunFile(runTask, count, 1);
            std::ofstream fw(sjmFile);
            fw << task;
            fw.close();
        }
        ++count;
        fr2.close();
    }
    fr1.close();
}
