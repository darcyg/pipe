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

void GenJob::genSpliterJob(const std::string& conf, Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/spliter";
    j->cmd.second += " -b " + mOpt->ioOpt.db_dir + "/barcode/barcode.conf";
    j->cmd.second += " -s " + conf;
    j->cmd.second += " -i " + lib1;
    j->cmd.second += " -I " + lib2;
    j->cmd.second += " -l " + j->workdir.second + "/" + j->pre + ".spliter.json";
    j->cmd.second += " -o " + j->workdir.second;
    j->memory.second = "5g";
    j->slots.second = "5";
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

void GenJob::genFqtoolJob(Job* j){
    std::string ofq1 = j->workdir.second + j->pre + ".R1.fq.gz";
    std::string ofq2 = j->workdir.second + j->pre + ".R2.fq.gz";
    j->cmd.second += mOpt->ioOpt.bin_dir + "/fqtool";
    j->cmd.second += " -q -a --detect_pe_adapter ";
    j->cmd.second += " -u --umi_location 6 --umi_skip_length 1 --umi_length 8 --umi_drop_comment";
    j->cmd.second += " -i " + lib1 + " -I " + lib2;
    j->cmd.second += " -o " + ofq1;
    j->cmd.second += " -O " + ofq2;
    j->cmd.second += " -J " + j->workdir.second + j->pre + ".fqtool.json";
    j->cmd.second += " -H " + j->workdir.second + j->pre + ".fqtool.html";
    j->memory.second = "1g";
    j->slots.second = "4";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";};
    j->o1 = ofq1;
    j->o2 = ofq2;
}

void GenJob::genSeqtkJob(Job* j){
    j->cmd.second = "rm -f " + j->workdir.second + util::basename(lib1);
    j->cmd.second += " && rm -f " + j->workdir.second + util::basename(lib2);
    j->cmd.second += " && ln -sf " + lib1 + " " + j->workdir.second;
    j->cmd.second += " && ln -sf " + lib2 + " " + j->workdir.second;
    if(mOpt->clOpt.dfq_vol != "VOL"){
        size_t vol = std::atoi(mOpt->clOpt.dfq_vol.c_str());
        if(vol == 0){
            vol = mOpt->minFqVolMap[j->optPre];
        }
        if(vol != 0){
            j->cmd.second = "rm -f " + j->workdir.second + util::basename(lib1);
            j->cmd.second += " && rm -f " + j->workdir.second + util::basename(lib2);
            j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/seqtk sample -2 -s 100 " + lib1 + " " + std::to_string(vol);
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

void GenJob::genFilterJob(Job* j){
    j->cmd.second += mOpt->ioOpt.bin_dir + "/filter -d";
    j->cmd.second += " -i " + lib1;
    j->cmd.second += " -I " + lib2;
    j->cmd.second += " -r " + mOpt->ioOpt.db_dir + "/ncrna/Homo_sapiens.GRCh37.ncrna.class.fa";
    j->cmd.second += " -o " + j->workdir.second + "/" + util::basename(lib1);
    j->cmd.second += " -O " + j->workdir.second + "/" + util::basename(lib2);
    j->cmd.second += " -l " + j->workdir.second + "/" + j->pre + ".filter.json";
    j->memory.second = "1g";
    j->slots.second = "6";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
    j->o1 = j->workdir.second + "/" + j->pre + ".R1.fq";
    j->o2 = j->workdir.second + "/" + j->pre + ".R2.fq";
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
    j->cmd.second += mOpt->ioOpt.bin_dir + "/duplexer";
    j->cmd.second += " -i " + bam;
    j->cmd.second += " -o " + j->workdir.second + "/" + j->pre + ".mkdup.bam";
    j->cmd.second += " -l " + j->workdir.second + "/" + j->pre + ".mkdup.json";
    j->cmd.second += " -r " + mOpt->clOpt.ref;
    j->cmd.second += " -m ";
    j->cmd.second += " > " + j->workdir.second + j->pre + ".mkdup.log 2>&1";
    j->cmd.second += " && " + mOpt->ioOpt.bin_dir + "/samtools sort -@ 8 -o " + j->workdir.second + j->pre + ".mkdup.sort.bam";
    j->cmd.second += " " + j->workdir.second + j->pre + ".mkdup.bam";
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
    j->cmd.second += " -b " + bam;
    j->cmd.second += " -r " + mOpt->clOpt.reg;
    j->cmd.second += " -o " + j->workdir.second + "/" + j->pre + ".bamqc.json";
    j->cmd.second += " -p " + j->pre;
    j->cmd.second += " > " + j->workdir.second + j->pre + ".log 2>&1";
    j->memory.second = "1g";
    j->slots.second = "4";
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
    std::string fqqclog = mOpt->ioOpt.cut_dir + "/" + j->pre + ".fqtool.json";
    std::string filtlog = mOpt->ioOpt.fil_dir + "/" + j->pre + ".filter.json";
    std::string bamqc = mOpt->ioOpt.bqc_dir + "/" + j->pre + ".bamqc.json";
    std::string fusrpt = mOpt->ioOpt.fus_dir + "/" + j->pre + "/" + j->pre + ".FusionReport.txt";
    std::string abundance = mOpt->ioOpt.exp_dir + "/" + j->pre + "/abundance.tsv";
    std::string ens2gen = mOpt->ioOpt.db_dir + "/NCBI/ensebml2genename";
    std::string outf = mOpt->ioOpt.rep_dir + "/" + j->pre + ".report.xlsx";
    j->cmd.second = mOpt->ioOpt.bin_dir + "/genrpt";
    j->cmd.second += " -q " + fqqclog + " -f " + filtlog + " -b " + bamqc + " -s " + fusrpt;
    j->cmd.second += " -k " + abundance + " -e " + ens2gen + " -g " + mOpt->clOpt.gset  + " -o " + outf;
    j->memory.second = "1g";
    j->slots.second = "1";
    j->sopt.second.append(" -l p=" + j->slots.second);
    j->sopt.second.append(" -l vf=" + j->memory.second);
    if(mOpt->clOpt.local){j->host.second = "localhost";}
}
