#include "genjob.h"

GenJob::GenJob(Options* opt){
    mOpt = opt;
}

GenJob::~GenJob(){
}

void GenJob::setLib(const std::string& l1, const std::string& l2){
    lib1 = l1;
    lib2 = l2;
}

void GenJob::genFastpJob(Job* j){
    std::string ofq1 = j->workdir.second + j->pre + ".R1.fq.gz";
    std::string ofq2 = j->workdir.second + j->pre + ".R2.fq.gz";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/fastp";
    j->cmd.second += " -i " + lib1 + " -I " + lib2;
    j->cmd.second += " -o " + ofq1;
    j->cmd.second += " -O " + ofq2;
    j->cmd.second += " -j " + j->workdir.second + j->pre + ".json";
    j->cmd.second += " -h " + j->workdir.second + j->pre + ".html";
    j->memory.second = "1g";
    j->slots.second = "2";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = ofq1;
    j->o2 = ofq2;
}
void GenJob::genSeqtkJob(Job* j){}
void GenJob::genFiltdbJob(Job* j){}
void GenJob::genAlignJob(Job* j){}
void GenJob::genMkdupJob(const std::string& bam, Job* j){}
void GenJob::genBamqcJob(const std::string& bam, Job* j){}
void GenJob::genFusionJob(Job* j){}
void GenJob::genExpressJob(Job* j){}
void GenJob::genReportJob(Job* j){}
