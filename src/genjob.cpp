#include "genjob.h"
#include "util.h"

GenJob::GenJob(Options* opt){
    mOpt = opt;
}

GenJob::~GenJob(){
}

void GenJob::setLib(const std::string& l1, const std::string& l2){
    lib1 = l1;
    lib2 = l2;
}

void GenJob::setBam(const std::string& b){
    bam = b;
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

void GenJob::genSplitrJob(const std::string& conf, Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/splitr";
    j->cmd.second += " -b " + mOpt->ioOpt.db_dir + "/barcode/barcode.conf";
    j->cmd.second += " -s " + conf;
    j->cmd.second += " -r " + lib1;
    j->cmd.second += " -R " + lib2;
    j->cmd.second += " -o " + j->workdir.second;
    j->cmd.second += " -p " + j->pre;
    j->memory.second = "4g";
    int count = 0;
    std::ifstream fr(conf.c_str());
    std::string tmpstr;
    while(std::getline(fr, tmpstr)){
        ++count;
    }
    fr.close();
    j->slots.second = std::to_string(count + 2);
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
}

void GenJob::genSeqtkJob(Job* j){
    j->cmd.second = "rm -f " + j->workdir.second + util::basename(lib1);
    j->cmd.second += " && rm -f " + j->workdir.second + util::basename(lib2);
    j->cmd.second += " && ln -sf " + lib1 + " " + j->workdir.second;
    j->cmd.second += " && ln -sf " + lib2 + " " + j->workdir.second;
    if(mOpt->clOpt.dfq_vol != "VOL"){
        size_t vol = std::atoi(mOpt->clOpt.dfq_vol.c_str());
        if(vol == 0){
            vol = getMinFqVol();
        }
        if(vol != 0){
            j->cmd.second = mOpt->ioOpt.bin_dir + "/seqtk sample -2 -s 100 " + lib1 + " " + std::to_string(vol);
            j->cmd.second += " > " + j->workdir.second + util::basename(lib1);
            j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/seqtk sample -2 -s 100 " + lib2 + " " + std::to_string(vol);
            j->cmd.second += " > " + j->workdir.second + util::basename(lib2);
        }
    }
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = j->workdir.second + util::basename(lib1);
    j->o2 = j->workdir.second + util::basename(lib2);
}

void GenJob::genFiltdbJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/filtdb";
    j->cmd.second += " -i " + lib1;
    j->cmd.second += " -I " + lib2;
    j->cmd.second += " -r " + mOpt->ioOpt.db_dir + "/rrna/human.rrna.fa";
    j->cmd.second += " -q -1 ";
    j->cmd.second += " -o " + j->workdir.second;
    j->cmd.second += " -p " + j->pre;
    j->memory.second = "1g";
    j->slots.second = "6";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
    j->o1 = j->workdir.second + util::basename(lib1);
    j->o2 = j->workdir.second + util::basename(lib2);
}

void GenJob::genAlignJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/bwa mem";
    j->cmd.second += " -C -t 8 " + mOpt->clOpt.ref + " " + lib1 + " " + lib2;
    j->cmd.second += " 2> " + j->workdir.second + j->pre + ".memalign.log";
    j->cmd.second += " | " + mOpt->ioOpt.bin_dir + "/samtools view -@ 8";
    j->cmd.second += " -o " + j->workdir.second + j->pre + ".aln.bam";
    j->cmd.second += " >> " + j->workdir.second + j->pre + ".memalign.log";
    j->cmd.second += " && rm -rf " + j->workdir.second + j->pre + ".aln.sort.bam*";
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools sort -@ 8";
    j->cmd.second += " -o " + j->workdir.second + j->pre + ".aln.sort.bam";
    j->cmd.second += " " + j->workdir.second + j->pre + ".aln.bam";
    j->cmd.second += " > " + j->workdir.second + j->pre + ".aln.sort.log 2>&1";
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools index ";
    j->cmd.second += j->workdir.second + j->pre + ".aln.sort.bam";
    j->cmd.second += " > " + j->workdir.second + j->pre + ".aln.sort.index.log 2>&1";
    j->memory.second = "16g";
    j->slots.second = "8";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = j->workdir.second + j->pre + ".aln.sort.bam";
}

void GenJob::genMkdupJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/mkdup";
    j->cmd.second += " -i " + bam;
    j->cmd.second += " -o " + j->workdir.second + util::basename(bam);
    j->cmd.second += " > " + j->workdir.second + j->pre + ".mkdup.log 2>&1";
    j->cmd.second += " && rm -rf " + j->workdir.second + j->pre + ".mkdup.sort.bam*";
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools sort -@ 8 -o " + j->workdir.second + j->pre + ".mkdup.sort.bam";
    j->cmd.second += " " + j->workdir.second + util::basename(bam);
    j->cmd.second += " > " + j->workdir.second + j->pre + ".mkdup.sort.log 2>&1";
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools index ";
    j->cmd.second += j->workdir.second + j->pre + ".mkdup.sort.bam";
    j->cmd.second += " > " + j->workdir.second + j->pre + ".mkdup.sort.index.log 2>&1";
    j->memory.second = "8g";
    j->slots.second = "8";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = j->workdir.second + j->pre + ".mkdup.sort.bam";
}

void GenJob::genBamqcJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/bamqc";
    j->cmd.second += " -i " + bam;
    j->cmd.second += " -b " + mOpt->clOpt.reg;
    j->cmd.second += " -o " + j->workdir.second;
    j->cmd.second += " -p " + j->pre;
    j->cmd.second += " -n " + j->pre;
    j->cmd.second += " > " + j->workdir.second + j->pre + "qc.log 2>&1";
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
}

void GenJob::genFusionJob(Job* j){
    j->cmd.second += "/share/public/software/mono/mono/bin/mono /share/public/software/oshell/oshell.exe";
    j->cmd.second += " --semap " + mOpt->ioOpt.db_dir + "/FusionMap/ Human.B37.3 RefGene";
    j->cmd.second += " " + j->workdir.second + j->pre + ".fs.cfg";
    j->cmd.second += " > " + j->workdir.second + j->pre + ".fs.log 2>&1";
    j->memory.second = "8g";
    j->slots.second = "8";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    FusionMapTemplate fm(lib1, lib2, j->workdir.second, j->pre);
    std::ofstream fw(j->workdir.second + j->pre + ".fs.cfg");
    fw << fm;
    fw.close();
}

void GenJob::genExpressJob(Job* j){
    j->cmd.second += "mkdir -p " + j->workdir.second + j->pre + " && ";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/kallisto";
    j->cmd.second += " quant -t 8 -i " + mOpt->ioOpt.db_dir + "/ensembl/Homo_sapiens.GRCh37.cdna.all.fa.idx";
    j->cmd.second += " -o " + j->workdir.second + j->pre + " ";
    j->cmd.second += lib1 + " " + lib2;
    j->memory.second = "4g";
    j->slots.second = "8";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
}

void GenJob::genReportJob(Job* j){
    std::string filtlog = mOpt->ioOpt.fil_dir + "/" + j->pre + ".filt.log";
    std::string ddpqc = mOpt->ioOpt.bqc_dir + "/" + j->pre + "_DDP_QC.tsv";
    std::string idpqc = mOpt->ioOpt.bqc_dir + "/" + j->pre + "_IDP_QC.tsv";
    std::string fusrpt = mOpt->ioOpt.fus_dir + "/" + j->pre + "/" + j->pre + ".FusionReport.txt";
    std::string abundance = mOpt->ioOpt.exp_dir + "/" + j->pre + "/abundance.tsv";
    std::string ens2gen = mOpt->ioOpt.db_dir + "/NCBI/ensebml2genename";
    std::string outf = mOpt->ioOpt.rep_dir + "/" + j->pre + ".report.xlsx";
    j->cmd.second = mOpt->ioOpt.bin_dir + "/genrpt";
    j->cmd.second += " -f " + filtlog + " -d " + ddpqc + " -i " + idpqc + " -s " + fusrpt;
    j->cmd.second += " -k " + abundance + " -e " + ens2gen + " -g " + mOpt->clOpt.gset  + " -o " + outf;
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
}

size_t GenJob::getMinFqVol(){
    std::vector<std::string> fname;
    util::list_dir(mOpt->ioOpt.spl_dir, fname);
    std::string logfile;
    for(size_t i = 0; i < fname.size(); ++i){
        if(util::ends_with(fname[i], ".spl.log")){
            logfile = util::joinpath(mOpt->ioOpt.spl_dir, fname[i]);
            break;
        }
    }
    if(logfile.empty()){
        return 0;
    }
    std::ifstream fr(logfile);
    std::string line;
    std::getline(fr, line);
    std::vector<size_t> readsNum;
    std::vector<std::string> vstr;
    while(std::getline(fr,line)){
        if(!util::starts_with(line, "Total/Get")){
            util::split(line, vstr,"\t");
            readsNum.push_back(std::atoi(vstr[1].c_str()));
        }
    }
    size_t minReadsNum = readsNum[0];
    for(size_t i = 1; i < readsNum.size(); ++i){
        minReadsNum = std::min(minReadsNum, readsNum[i]);
    }
    return minReadsNum;
}
