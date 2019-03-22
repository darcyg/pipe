#include "sjm.h"
#include "dirutil.h"

// show stage_marker<->ana_stage relationship
void sjm::show_mark(){
    std::vector<std::string> stg = {"cutadaptor by fastp", "split read by splitr", "filter rrna by filtdb",
                                    "downsample fastq by seqtk", "genome slignment by bwa", "markdup by mkdp",
                                    "bam QC by bamqc", "fusion calling by fusionMap", "express quant by kallisto",
                                    "report by genrep"};
    std::cout << std::left;
    std::cout << "  ┌--------------------------------------┐" << std::endl;
    std::cout << "  |" << std::setw(10) << "Marker" << std::setw(28) << "Analysis" << "|" << std::endl;
    for(size_t i = 0; i < stg.size(); ++i){
        std::cout << "  |--------------------------------------|" << std::endl;
        std::cout << "  |" << std::setw(10) << i + 1 << std::setw(28) << stg[i] << "|" << std::endl;
    }
    std::cout << "  └--------------------------------------┘" << std::endl;
}

// update args after commandline args parsed
void sjm::update_args(sjm::args& a){
    a.bin_dir = util::dirname(dirutil::getExecutablePath()); 
    a.sample_list = util::abspath(a.sample_list);
    a.out_dir = util::abspath(a.out_dir);
    a.db_dir = util::dirname(a.bin_dir) + "/db/";
    if(a.ana_marker.empty()){
        for(int i = a.ini_marker; i <= a.end_marker; ++i){
            a.ana_marker.push_back(i);
        }
    }else{
        std::remove_if(a.ana_marker.begin(), a.ana_marker.end(), [&a](int& e){return e < a.ini_marker || e > a.end_marker;});
    }
    if(a.gset.empty()){
        a.gset = a.db_dir + "gset/kgset.tsv";
    }
    std::string sep = "/";
    a.sjm_dir = a.out_dir + sep + a.sjm_dir;
    a.cut_dir = a.out_dir + sep + a.cut_dir;
    a.spl_dir = a.out_dir + sep + a.spl_dir;
    a.fil_dir = a.out_dir + sep + a.fil_dir;
    a.dfq_dir = a.out_dir + sep + a.dfq_dir;
    a.aln_dir = a.out_dir + sep + a.aln_dir;
    a.mkd_dir = a.out_dir + sep + a.mkd_dir;
    a.bqc_dir = a.out_dir + sep + a.bqc_dir;
    a.fus_dir = a.out_dir + sep + a.fus_dir;
    a.exp_dir = a.out_dir + sep + a.exp_dir;
    a.rep_dir = a.out_dir + sep + a.rep_dir;
    a.log_dir = a.out_dir + sep + a.log_dir;
}

// prepare output sub directories
void sjm::gen_dir(const sjm::args& a){
    util::makedir(a.sjm_dir);
    util::makedir(a.cut_dir);
    util::makedir(a.spl_dir);
    util::makedir(a.fil_dir);
    util::makedir(a.dfq_dir);
    util::makedir(a.aln_dir);
    util::makedir(a.mkd_dir);
    util::makedir(a.bqc_dir);
    util::makedir(a.fus_dir);
    util::makedir(a.exp_dir);
    util::makedir(a.rep_dir);
    util::makedir(a.log_dir);
}

// generate fastp job
void sjm::gen_fastp_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    std::string ofq1 = j.workdir.second + pre + ".R1.fq.gz";
    std::string ofq2 = j.workdir.second + pre + ".R2.fq.gz";
    j.cmd.second += a.bin_dir + "/fastp";
    j.cmd.second += " -i " + lib1 + " -I " + lib2; 
    j.cmd.second += " -o " + ofq1;
    j.cmd.second += " -O " + ofq2;
    j.cmd.second += " -j " + j.workdir.second + pre + ".json";
    j.cmd.second += " -h " + j.workdir.second + pre + ".html";
    j.memory.second = "1g";
    j.slots.second = "2";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";};
    j.o1 = ofq1;
    j.o2 = ofq2;
}

// generate splitr job
void sjm::gen_splitr_job(const sjm::args& a, const std::string& lib1, const std::string& lib2,  const std::string& sconf, const std::string& pre, sjm::job& j){
    j.cmd.second += a.bin_dir + "/splitr";
    j.cmd.second += " -b " + a.db_dir + "/barcode/barcode.conf";
    j.cmd.second += " -s " + sconf;
    j.cmd.second += " -r " + lib1;
    j.cmd.second += " -R " + lib2;
    j.cmd.second += " -o " + j.workdir.second;
    j.cmd.second += " -p " + pre;
    j.memory.second = "4g";
    int count = 0;
    std::ifstream fr(sconf.c_str());
    std::string tmpstr;
    while(std::getline(fr, tmpstr)){
        ++count;
    }
    fr.close();
    j.slots.second = std::to_string(count + 2);
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";};
}

// generate filtdb job
void sjm::gen_filtdb_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.cmd.second += a.bin_dir + "/filtdb";
    j.cmd.second += " -i " + lib1;
    j.cmd.second += " -I " + lib2;
    j.cmd.second += " -r " + a.db_dir + "/rrna/human.rrna.fa";
    j.cmd.second += " -q -1 ";
    j.cmd.second += " -o " + j.workdir.second;
    j.cmd.second += " -p " + pre;
    j.memory.second = "1g";
    j.slots.second = "6";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";}
    j.o1 = j.workdir.second + util::basename(lib1);
    j.o2 = j.workdir.second + util::basename(lib2);
}

// generate sektk job
void sjm::gen_seqtk_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    if(a.dfq_vol == "VOL"){
        j.cmd.second = "rm -f " + j.workdir.second + util::basename(lib1);
        j.cmd.second += " && rm -f " + j.workdir.second + util::basename(lib2);
        j.cmd.second += " && ln -sf " + lib1 + " " + j.workdir.second;
        j.cmd.second += " && ln -sf " + lib2 + " " + j.workdir.second;
    }
    else{
        size_t vol = std::atoi(a.dfq_vol.c_str());
        if(vol == 0){
            vol = sjm::get_min_sublib_vol(a);
        }
        j.cmd.second = a.bin_dir + "/seqtk sample -2 -s 100 " + lib1 + " " + std::to_string(vol);
        j.cmd.second += " > " + j.workdir.second + util::basename(lib1);
        j.cmd.second += " && " + a.bin_dir + "/seqtk sample -2 -s 100 " + lib2 + " " + std::to_string(vol);
        j.cmd.second += " > " + j.workdir.second + util::basename(lib2);
    }
    j.memory.second = "1g";
    j.slots.second = "1";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";};
    j.o1 = j.workdir.second + util::basename(lib1);
    j.o2 = j.workdir.second + util::basename(lib2);
}

// generate alignment job
void sjm::gen_align_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.cmd.second += a.bin_dir + "/bwa mem";
    j.cmd.second += " -C -t 8 " + a.ref + " " + lib1 + " " + lib2;
    j.cmd.second += " 2> " + j.workdir.second + pre + ".memalign.log";
    j.cmd.second += " | " + a.bin_dir + "/samtools view -@ 8";
    j.cmd.second += " -o " + j.workdir.second + pre + ".aln.bam";
    j.cmd.second += " >> " + j.workdir.second + pre + ".memalign.log";
    j.cmd.second += " && rm -rf " + j.workdir.second + pre + ".aln.sort.bam*";
    j.cmd.second += " && " + a.bin_dir + "/samtools sort -@ 8";
    j.cmd.second += " -o " + j.workdir.second + pre + ".aln.sort.bam";
    j.cmd.second += " " + j.workdir.second + pre + ".aln.bam";
    j.cmd.second += " > " + j.workdir.second + pre + ".aln.sort.log 2>&1";
    j.cmd.second += " && " + a.bin_dir + "/samtools index ";
    j.cmd.second += j.workdir.second + pre + ".aln.sort.bam";
    j.cmd.second += " > " + j.workdir.second + pre + ".aln.sort.index.log 2>&1";
    j.memory.second = "16g";
    j.slots.second = "8";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";};
    j.o1 = j.workdir.second + pre + ".aln.sort.bam";
}

// generate mkdup job
void sjm::gen_mkdup_job(const sjm::args& a, const std::string& bam, const std::string& pre, sjm::job& j){
    j.cmd.second += a.bin_dir + "/mkdup";
    j.cmd.second += " -i " + bam;
    j.cmd.second += " -o " + j.workdir.second + util::basename(bam);
    j.cmd.second += " > " + j.workdir.second + pre + ".mkdup.log 2>&1";
    j.cmd.second += " && rm -rf " + j.workdir.second + pre + ".mkdup.sort.bam*";
    j.cmd.second += " && " + a.bin_dir + "/samtools sort -@ 8 -o " + j.workdir.second + pre + ".mkdup.sort.bam";
    j.cmd.second += " " + j.workdir.second + util::basename(bam);
    j.cmd.second += " > " + j.workdir.second + pre + ".mkdup.sort.log 2>&1";
    j.cmd.second += " && " + a.bin_dir + "/samtools index ";
    j.cmd.second += j.workdir.second + pre + ".mkdup.sort.bam";
    j.cmd.second += " > " + j.workdir.second + pre + ".mkdup.sort.index.log 2>&1";
    j.memory.second = "8g";
    j.slots.second = "8";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";};
    j.o1 = j.workdir.second + pre + ".mkdup.sort.bam";
}

// generate bamqc job
void sjm::gen_bamqc_job(const sjm::args& a, const std::string& bam, const std::string& pre, sjm::job& j){
    j.cmd.second += a.bin_dir + "/bamqc";
    j.cmd.second += " -i " + bam;
    j.cmd.second += " -b " + a.reg;
    j.cmd.second += " -o " + j.workdir.second;
    j.cmd.second += " -p " + pre;
    j.cmd.second += " -n " + pre;
    j.cmd.second += " > " + j.workdir.second + pre + "qc.log 2>&1";
    j.memory.second = "1g";
    j.slots.second = "1";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";};
}

// generate fusion job
void sjm::gen_fusion_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.cmd.second += "/share/public/software/mono/mono/bin/mono /share/public/software/oshell/oshell.exe";
    j.cmd.second += " --semap " + a.db_dir + "/FusionMap/ Human.B37.3 RefGene";
    j.cmd.second += " " + j.workdir.second + pre + ".fs.cfg";
    j.cmd.second += " > " + j.workdir.second + pre + ".fs.log 2>&1";
    j.memory.second = "8g";
    j.slots.second = "8";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";};
    sjm::fmap fm(lib1, lib2, j.workdir.second, pre);
    std::ofstream fw(j.workdir.second + pre + ".fs.cfg");
    fw << fm;
    fw.close();
}

// generate express job
void sjm::gen_express_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.cmd.second += "mkdir -p " + j.workdir.second + pre + " && ";
    j.cmd.second += a.bin_dir + "/kallisto";
    j.cmd.second += " quant -t 8 -i " + a.db_dir + "/ensembl/Homo_sapiens.GRCh37.cdna.all.fa.idx";
    j.cmd.second += " -o " + j.workdir.second + pre + " ";
    j.cmd.second += lib1 + " " + lib2;
    j.memory.second = "4g";
    j.slots.second = "8";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";}
}

// generate report job
void sjm::gen_report_job(const sjm::args& a, const std::string& lib, sjm::job& j){
    std::string filtlog = a.fil_dir + "/" + lib + ".filt.log";
    std::string ddpqc = a.bqc_dir + "/" + lib + "_DDP_QC.tsv";
    std::string idpqc = a.bqc_dir + "/" + lib + "_IDP_QC.tsv";
    std::string fusrpt = a.fus_dir + "/" + lib + "/" + lib + ".FusionReport.txt";
    std::string abundance = a.exp_dir + "/" + lib + "/abundance.tsv";
    std::string ens2gen = a.db_dir + "/NCBI/ensebml2genename";
    std::string outf = a.rep_dir + "/" + lib + ".report.xlsx";
    j.cmd.second = a.bin_dir + "/genrpt";
    j.cmd.second += " -f " + filtlog + " -d " + ddpqc + " -i " + idpqc + " -s " + fusrpt;
    j.cmd.second += " -k " + abundance + " -e " + ens2gen + " -g " + a.gset  + " -o " + outf;
    j.memory.second = "1g";
    j.slots.second = "1";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    if(a.local){j.host.second = "localhost";}
}

// generate prelib task
void sjm::gen_prelib_task(const sjm::args& a, sjm::pipeline& p){
    std::ifstream fr(a.sample_list);
    std::string sample_no, flow_cell, lib_name, read1, read2, barcode_conf, tmp_str;
    std::getline(fr, tmp_str);
    std::istringstream iss;
    int count = 0;
    while(std::getline(fr, tmp_str)){
        iss.clear();
        iss.str(tmp_str);
        iss >> sample_no >> flow_cell >> lib_name >> read1 >> read2 >> barcode_conf;
        sjm::job jfastp("fastp", a.cut_dir, lib_name, 1);
        sjm::job jsplitr("splitr", a.spl_dir, lib_name, 2);
        sjm::gen_fastp_job(a, read1, read2, lib_name, jfastp);
        sjm::gen_splitr_job(a, jfastp.o1, jfastp.o2, barcode_conf, lib_name, jsplitr);
        sjm::task t;
        t.joblist.resize(2);
        t.joblist[0].push_back(jfastp);
        t.joblist[1].push_back(jsplitr);
        std::string sjmstatus(a.sjm_dir + "/" + lib_name + "_prelib.sjm.status");
        std::map<std::string, std::string> jmap;
        sjm::get_status(jmap, sjmstatus);
        for(auto& e: t.joblist){
            for(auto& f: e){
                if(std::find(a.ana_marker.cbegin(), a.ana_marker.cend(), f.stage_marker) == a.ana_marker.cend()){
                    f.status.second = "done";
                }
                if(a.update && jmap.find(f.name.second) != jmap.end()){
                    f.status.second = jmap[f.name.second];
                }
                if(!a.queue.empty()){
                    f.queue.second = a.queue;
                }
            }
        }
        std::string sjmfile(a.sjm_dir + "/" + lib_name + "_prelib.sjm");
        std::string sjmlogf(a.sjm_dir + "/" + lib_name + "_prelib.log");
        sjm::runfile rf;
        rf.sjmcmd = a.bin_dir + "/sjm -i --log_level verbose -l " + sjmlogf + " " + sjmfile;
        rf.sucf = a.log_dir + "/" + lib_name + ".Prelib.SUCCESS";
        rf.faif = a.log_dir + "/" + lib_name + ".Prelib.FAIL";
        rf.logf = sjmlogf;
        p.pipelist[count][0].push_back(rf);
        ++count;
        std::ofstream fw(sjmfile);
        fw << t;
        fw.close();
    }
    fr.close();
}

// generate analib task
void sjm::gen_analib_task(const sjm::args& a, sjm::pipeline& p){
    std::ifstream fr1(a.sample_list);
    std::string sample_no, flow_cell, lib_name, read1, read2, barcode_conf;
    std::string sublib, subr1, subr2, tmp_str;
    std::getline(fr1, tmp_str);
    std::istringstream iss1, iss2;
    int count = 0;
    while(std::getline(fr1, tmp_str)){
        iss1.clear();
        iss1.str(tmp_str);
        iss1 >> sample_no >> flow_cell >> lib_name >> read1 >> read2 >> barcode_conf;
        std::ifstream fr2(barcode_conf);
        while(std::getline(fr2, tmp_str)){
            iss2.clear();
            iss2.str(tmp_str);
            iss2 >> sublib;
            subr1 = a.spl_dir + "/" + sublib + ".R1.fq";
            subr2 = a.spl_dir + "/" + sublib + ".R2.fq";
            sjm::job jfiltdb("filtdb", a.fil_dir, sublib, 3);
            sjm::job jseqtk("seqtk", a.dfq_dir, sublib, 4);
            sjm::job jaln("bwa", a.aln_dir, sublib, 5);
            sjm::job jmkdup("mkdup", a.mkd_dir, sublib, 6);
            sjm::job jbamqc("bamqc", a.bqc_dir, sublib, 7);
            sjm::job jfusion("fusionmap", a.fus_dir, sublib, 8);
            sjm::job jexpress("kalliso", a.exp_dir, sublib, 9);
            sjm::job jreport("genrpt", a.rep_dir, sublib, 10);
            sjm::gen_filtdb_job(a, subr1, subr2, sublib, jfiltdb);
            sjm::gen_seqtk_job(a, jfiltdb.o1, jfiltdb.o2, sublib, jseqtk);
            sjm::gen_align_job(a, jseqtk.o1, jseqtk.o2, sublib, jaln);
            sjm::gen_mkdup_job(a, jaln.o1, sublib, jmkdup);
            sjm::gen_bamqc_job(a, jmkdup.o1, sublib, jbamqc);
            sjm::gen_fusion_job(a, jseqtk.o1, jseqtk.o2, sublib, jfusion);
            sjm::gen_express_job(a, jseqtk.o1, jseqtk.o2, sublib, jexpress);
            sjm::gen_report_job(a, sublib, jreport);
            sjm::task t;
            t.joblist.resize(6);
            t.joblist[0].push_back(jfiltdb);
            t.joblist[1].push_back(jseqtk);
            t.joblist[2].push_back(jaln);
            t.joblist[2].push_back(jfusion);
            t.joblist[2].push_back(jexpress);
            t.joblist[3].push_back(jmkdup);
            t.joblist[4].push_back(jbamqc);
            t.joblist[5].push_back(jreport);
            std::string sjmstatus = a.sjm_dir + "/" + sublib + "_analib.sjm.status";
            std::map<std::string, std::string> jmap;
            sjm::get_status(jmap, sjmstatus);
            for(auto& e: t.joblist){
                for(auto& f: e){
                    if(std::find(a.ana_marker.cbegin(), a.ana_marker.cend(), f.stage_marker) == a.ana_marker.cend()){
                        f.status.second = "done";
                    }
                    if(a.update && jmap.find(f.name.second) != jmap.end()){
                        f.status.second = jmap[f.name.second];
                    }
                    if(!a.queue.empty()){
                        f.queue.second = a.queue;
                    }
                }
            }
            std::string sjmfile(a.sjm_dir + "/" + sublib + "_analib.sjm");
            std::string sjmlogf(a.sjm_dir + "/" + sublib + "_analib.log");
            sjm::runfile rf;
            rf.sjmcmd = a.bin_dir + "/sjm -i --log_level verbose -l " + sjmlogf + " " + sjmfile;
            rf.sucf = a.log_dir + "/" + sublib + ".Analysis.SUCCESS";
            rf.faif = a.log_dir + "/" + sublib + ".Analysis.FAIL";
            rf.logf = sjmlogf;
            p.pipelist[count][1].push_back(rf);
            std::ofstream fw(sjmfile);
            fw << t;
            fw.close();
        }
        ++count;
        fr2.close();
    }
    fr1.close();
}


// initialize pipeline
void sjm::init_pipeline(const sjm::args& a, sjm::pipeline& p){
    std::ifstream fr(a.sample_list);
    std::string tmp_str;
    int c = 0;
    while(std::getline(fr, tmp_str)){
        ++c;
    }
    // c - 1 parallel pipeline
    p.pipelist.resize(c - 1);
    for(auto& e: p.pipelist){
        // each pipeline has two tasks
        e.resize(2);
    }
    p.faif = a.log_dir + "/FAIL";
    p.sucf = a.log_dir + "/SUCCESS";
    p.ret = 0;
}

// running only one sjm file
int sjm::run_sjm(const std::string& sjm){
    int ret = std::system(sjm.c_str());
    return ret;
}

// running only one parallel task
int sjm::run_task(std::vector<std::vector<sjm::runfile>>& t){
    std::vector<std::future<int>> fv;
    int fret = 0;
    for(size_t i = 0; i < t.size(); ++i){
        fv.clear();
        for(size_t j = 0; j < t[i].size(); ++j){
            fv.push_back(std::async(std::launch::async, sjm::run_sjm, std::ref(t[i][j].sjmcmd)));
        }
        for(size_t k = 0; k < t[i].size(); ++k){
            t[i][k].ret = fv[k].get();
            if(t[i][k].ret){
                remove(t[i][k].sucf.c_str());
                std::ofstream fw(t[i][k].faif.c_str());
                fw.close();
                fret = 1;
            }else{
                remove(t[i][k].faif.c_str());
                std::ofstream fw(t[i][k].sucf.c_str());
                fw.close();
            }
        }
    }
    return fret;
}

// running pipeline
int sjm::pipeline::run_pipe(){
    std::vector<std::future<int>> fv;
    for(size_t i = 0; i < this->pipelist.size(); ++i){
        fv.push_back(std::async(std::launch::async, sjm::run_task, std::ref(this->pipelist[i])));
    }
    for(auto& e: fv){
        this->ret += std::abs(e.get());
    }
    if(this->ret){
        remove(this->sucf.c_str());
        std::ofstream fw(this->faif.c_str());
        fw.close();
    }else{
        remove(this->faif.c_str());
        std::ofstream fw(this->sucf.c_str());
        fw.close();
    }
    return this->ret;
}

// prepare to rerun
void sjm::pipeline::pre_rerun(){
    std::string orisjm, newsjm, baksjm;
    for(size_t i = 0; i < this->pipelist.size(); ++i){
        for(size_t j = 0; j < this->pipelist[i].size(); ++j){
            for(auto& e :this->pipelist[i][j]){
                remove(e.faif.c_str());
                remove(e.sucf.c_str());
                remove(e.logf.c_str());
                auto p = e.sjmcmd.find_last_of(" "); //get the position just befor sjm path
                orisjm = e.sjmcmd.substr(p + 1);
                newsjm = orisjm + ".status";
                baksjm = newsjm + ".bak";
                remove(baksjm.c_str());
                std::ifstream fr(newsjm.c_str());
                if(fr.is_open()){
                    fr.close();
                    if(!this->update){
                        std::rename(newsjm.c_str(), orisjm.c_str());
                    }
                }
            }
        }
    }
}

// test logfile to see any error happen or not
bool sjm::test_job_fail(const std::string& log){
    std::ifstream fr(log);
    std::ostringstream oss;
    oss << fr.rdbuf();
    std::string logs = oss.str();
    auto p = logs.find("Failed jobs");
    fr.close();
    return p != std::string::npos;
}

// get map of stage finish status
void sjm::get_status(std::map<std::string, std::string>& jmap, std::string& jfile){
    std::pair<std::string, std::string> pjob;
    std::ifstream fr(jfile);
    std::istringstream iss;
    std::string line, prefix, suffix;
    while(std::getline(fr, line)){
        iss.clear();
        iss.str(line);
        iss >> prefix >> suffix;
        if(prefix == "name"){
            pjob.first = suffix;
        }
        if(prefix == "status"){
            pjob.second = suffix;
            jmap[pjob.first] = pjob.second;
        }
    }
    fr.close();
}

// get split sublibrary minimum reads
size_t sjm::get_min_sublib_vol(const sjm::args& a){
    tinydir_dir dir;
    tinydir_open(&dir, a.spl_dir.c_str());
    std::string logfile;
    while(dir.has_next){
        tinydir_file file;
        tinydir_readfile(&dir, &file);
        logfile = std::string(file.name);
        if(util::ends_with(logfile, ".spl.log")){
            break;
        }
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
    for(int i = 1; i < readsNum.size(); ++i){
        minReadsNum = std::min(minReadsNum, readsNum[i]);
    }
    return minReadsNum;
}
    
