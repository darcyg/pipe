#include "sjm.h"
#include "util.h"

// prepare output sub directories
void sjm::gen_dir(const sjm::args& a){
    std::string sep = "/";
    util::make_dirs(a.out_dir + sep + a.sjm_dir);
    util::make_dirs(a.out_dir + sep + a.cut_dir);
    util::make_dirs(a.out_dir + sep + a.spl_dir);
    util::make_dirs(a.out_dir + sep + a.fil_dir);
    util::make_dirs(a.out_dir + sep + a.dfq_dir);
    util::make_dirs(a.out_dir + sep + a.aln_dir);
    util::make_dirs(a.out_dir + sep + a.mkd_dir);
    util::make_dirs(a.out_dir + sep + a.bqc_dir);
    util::make_dirs(a.out_dir + sep + a.fus_dir);
    util::make_dirs(a.out_dir + sep + a.rep_dir);
    util::make_dirs(a.out_dir + sep + a.log_dir);
}

// generate fastp job
void sjm::gen_fastp_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.name.second = "fastp";
    j.workdir.second = a.out_dir + "/" + a.cut_dir + "/";
    std::string ofq1 = j.workdir.second + pre + ".R1.fq.gz";
    std::string ofq2 = j.workdir.second + pre + ".R2.fq.gz";
    std::string sgee = j.workdir.second + pre + ".sub.e";
    std::string sgeo = j.workdir.second + pre + ".sub.o";
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
    j.sopt.second.append(" -e " + sgee + " -o " + sgeo);
    if(a.local){j.host.second = "localhost";};
    j.o1 = ofq1;
    j.o2 = ofq2;
}

// generate splitr job
void sjm::gen_splitr_job(const sjm::args& a, const std::string& lib1, const std::string& lib2,  const std::string& sconf, const std::string& pre, sjm::job& j){
    j.name.second = "splitr";
    j.workdir.second = a.out_dir + "/" + a.spl_dir + "/";
    std::string sgee = j.workdir.second + pre + ".sub.e";
    std::string sgeo = j.workdir.second + pre + ".sub.o";
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
    j.sopt.second.append(" -e " + sgee + " -o " + sgeo);
    if(a.local){j.host.second = "localhost";};
}

// generate filtdb job
void sjm::gen_filtdb_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.name.second = "filtdb";
    j.workdir.second = a.out_dir + "/" + a.fil_dir + "/";
    std::string sgee = j.workdir.second + pre + ".sub.e";
    std::string sgeo = j.workdir.second + pre + ".sub.o";
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
    j.sopt.second.append(" -e " + sgee + " -o " + sgeo);
    if(a.local){j.host.second = "localhost";}
    j.o1 = j.workdir.second + util::get_basename(lib1);
    j.o2 = j.workdir.second + util::get_basename(lib2);
}

// generate sektk job
void sjm::gen_seqtk_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.name.second = "seqtk";
    j.workdir.second = a.out_dir + "/" + a.dfq_dir + "/";
    std::string sgee = j.workdir.second + pre + ".sub.e";
    std::string sgeo = j.workdir.second + pre + ".sub.o";
    if(a.dfq_vol == "VOL"){
        j.cmd.second = "ln -sf " + lib1 + " " + j.workdir.second;
        j.cmd.second += " && ln -sf " + lib2 + " " + j.workdir.second;
        j.cmd.second += " > " + j.workdir.second + pre + ".link.log 2>&1";
    }else{
        j.cmd.second = a.bin_dir + "/seqtk sample -2 -s 100 " + lib1 + " " + a.dfq_vol;
        j.cmd.second += " > " + j.workdir.second + util::get_basename(lib1);
        j.cmd.second += " && " + a.bin_dir + "/seqtk sample -2 -s 100 " + lib2 + " " + a.dfq_vol;
        j.cmd.second += " > " + j.workdir.second + util::get_basename(lib2);
    }
    j.memory.second = "1g";
    j.slots.second = "1";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    j.sopt.second.append(" -e " + sgee + " -o " + sgeo);
    if(a.local){j.host.second = "localhost";};
    j.o1 = j.workdir.second + util::get_basename(lib1);
    j.o2 = j.workdir.second + util::get_basename(lib2);
}

// generate alignment job
void sjm::gen_align_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.name.second = "bwa";
    j.workdir.second = a.out_dir + "/" + a.aln_dir + "/";
    std::string sgee = j.workdir.second + pre + ".sub.e";
    std::string sgeo = j.workdir.second + pre + ".sub.o";
    j.cmd.second += a.bin_dir + "/bwa mem";
    j.cmd.second += " -C -t 8 " + a.db_dir + "/refMrna/refMrna.fa " + lib1 + " " + lib2;
    j.cmd.second += " 2> " + j.workdir.second + pre + ".memalign.log";
    j.cmd.second += " | " + a.bin_dir + "/samtools view -@ 8";
    j.cmd.second += " -o " + j.workdir.second + pre + ".aln.bam";
    j.cmd.second += " >> " + j.workdir.second + pre + ".memalign.log";
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
    j.sopt.second.append(" -e " + sgee + " -o " + sgeo);
    if(a.local){j.host.second = "localhost";};
    j.o1 = j.workdir.second + pre + ".aln.sort.bam";
}

// generate mkdup job
void sjm::gen_mkdup_job(const sjm::args& a, const std::string& bam, const std::string& pre, sjm::job& j){
    j.name.second = "mkdup";
    j.workdir.second = a.out_dir + "/" + a.mkd_dir + "/";
    std::string sgee = j.workdir.second + pre + ".sub.e";
    std::string sgeo = j.workdir.second + pre + ".sub.o";
    j.cmd.second += a.bin_dir + "/mkdup";
    j.cmd.second += " -i " + bam;
    j.cmd.second += " -o " + j.workdir.second + util::get_basename(bam);
    j.cmd.second += " > " + j.workdir.second + pre + ".mkdup.log 2>&1";
    j.cmd.second += " && " + a.bin_dir + "/samtools sort -@ 8 -o " + j.workdir.second + pre + ".mkdup.sort.bam";
    j.cmd.second += " > " + j.workdir.second + pre + ".mkdup.sort.log 2>&1";
    j.cmd.second += " && " + a.bin_dir + "/samtools index ";
    j.cmd.second += j.workdir.second + pre + ".mkdup.sort.bam";
    j.cmd.second += " > " + j.workdir.second + pre + ".mkdup.sort.index.log 2>&1";
    j.memory.second = "8g";
    j.slots.second = "8";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    j.sopt.second.append(" -e " + sgee + " -o " + sgeo);
    if(a.local){j.host.second = "localhost";};
    j.o1 = j.workdir.second + pre + ".mkdup.sort.bam";
}

// generate bamqc job
void sjm::gen_bamqc_job(const sjm::args& a, const std::string& bam, const std::string& pre, sjm::job& j){
    j.name.second = "bamqc";
    j.workdir.second = a.out_dir + "/" + a.bqc_dir + "/";
    std::string sgee = j.workdir.second + pre + ".sub.e";
    std::string sgeo = j.workdir.second + pre + ".sub.o";
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
    j.sopt.second.append(" -e " + sgee + " -o " + sgeo);
    if(a.local){j.host.second = "localhost";};
}

// generate fusion job
void sjm::gen_fusion_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j){
    j.name.second = "fusion";
    j.workdir.second = a.out_dir + "/" + a.fus_dir + "/";
    std::string sgee = j.workdir.second + pre + ".sub.e";
    std::string sgeo = j.workdir.second + pre + ".sub.o";
    j.cmd.second += "/share/public/software/mono/mono/bin/mono /share/public/software/oshell/oshell.exe";
    j.cmd.second += " --semap " + a.db_dir + "/FusionMap/ Human.B37.3 RefGene";
    j.cmd.second += " " + j.workdir.second + pre + ".fs.cfg";
    j.cmd.second += " > " + j.workdir.second + pre + ".fs.log 2>&1";
    j.memory.second = "8g";
    j.slots.second = "8";
    j.sopt.second.append(" -l p=" + j.slots.second);
    j.sopt.second.append(" -l vf=" + j.memory.second);
    j.sopt.second.append(" -e " + sgee + " -o " + sgeo);
    if(a.local){j.host.second = "localhost";};
    sjm::fmap fm(lib1, lib2, j.workdir.second, pre);
    std::ofstream fw(j.workdir.second + pre + ".fs.cfg");
    fw << fm;
    fw.close();
}

// generate report job
void sjm::gen_report_job(const sjm::args& a, const std::string& lib){
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
        sjm::job jfastp, jsplitr;
        sjm::gen_fastp_job(a, read1, read2, lib_name, jfastp);
        sjm::gen_splitr_job(a, jfastp.o1, jfastp.o2, barcode_conf, lib_name, jsplitr);
        sjm::task t;
        t.joblist.resize(2);
        t.joblist[0].push_back(jfastp);
        t.joblist[1].push_back(jsplitr);
        std::string sjmfile(a.out_dir + "/" + a.sjm_dir + "/" + lib_name + "_prelib.sjm");
        p.pipelist[count][0].push_back(a.bin_dir + "/sjm -i " + sjmfile);
        ++p.njob[0];
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
            subr1 = a.out_dir + "/" + a.spl_dir + "/" + sublib + ".R1.fq";
            subr2 = a.out_dir + "/" + a.spl_dir + "/" + sublib + ".R2.fq";
            sjm::job jfiltdb, jseqtk, jaln, jmkdup, jbamqc, jfusion;
            sjm::gen_filtdb_job(a, subr1, subr2, sublib, jfiltdb);
            sjm::gen_seqtk_job(a, jfiltdb.o1, jfiltdb.o2, sublib, jseqtk);
            sjm::gen_align_job(a, jseqtk.o1, jseqtk.o2, sublib, jaln);
            sjm::gen_mkdup_job(a, jaln.o1, sublib, jmkdup);
            sjm::gen_bamqc_job(a, jmkdup.o1, sublib, jbamqc);
            sjm::gen_fusion_job(a, jseqtk.o1, jseqtk.o2, sublib, jfusion);
            sjm::task t;
            t.joblist.resize(6);
            t.joblist[0].push_back(jfiltdb);
            t.joblist[1].push_back(jseqtk);
            t.joblist[2].push_back(jaln);
            t.joblist[2].push_back(jfusion);
            t.joblist[3].push_back(jmkdup);
            t.joblist[4].push_back(jbamqc);
            std::string sjmfile(a.out_dir + "/" + a.sjm_dir + "/" + sublib + "_analib.sjm");
            p.pipelist[count][1].push_back(a.bin_dir + "/sjm -i " + sjmfile);
            ++p.njob[1];
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
    p.njob.resize(2);
    p.njob[0] = p.njob[1] = 0;
}

// running only one sjm file
int sjm::run_sjm(const std::string& sjm){
    int ret = std::system(sjm.c_str());
    return ret;
}

// running only one parallel task
int sjm::run_task(const std::vector<std::vector<std::string>>& t){
    std::vector<std::future<int>> fv;
    int ret = 0;
    int fret = 0;
    for(size_t i = 0; i < t.size(); ++i){
        fv.clear();
        for(size_t j = 0; j < t[i].size(); ++j){
            fv.push_back(std::async(std::launch::async, sjm::run_sjm, t[i][j]));
        }
        for(auto& e: fv){
            ret = e.get();
            if(ret != 0){
                fret = ret;
            }
        }
    }
    return fret;
}

// running pipeline
int sjm::pipeline::run_pipe(){
    int ret, fret = 0;
    std::vector<std::future<int>> fv;
    for(size_t i = 0; i < this->pipelist.size(); ++i){
        fv.push_back(std::async(std::launch::async, sjm::run_task, this->pipelist[i]));
    }
    for(auto& e: fv){
        ret = e.get();
        if(ret != 0){
            fret = ret;
        }
    }
    return fret;
}

// prepare to rerun
void sjm::pipeline::pre_rerun(){
    std::string orisjm, newsjm, baksjm;
    for(size_t i = 0; i < this->pipelist.size(); ++i){
        for(size_t j = 0; j < this->pipelist[i].size(); ++j){
            for(auto& e :this->pipelist[i][j]){
                auto p = e.find_last_of(" "); //get the position just befor sjm path
                orisjm = e.substr(p + 1);
                newsjm = orisjm + ".status";
                baksjm = newsjm + ".bak";
                remove(baksjm.c_str());
                std::ifstream fr(newsjm.c_str());
                if(fr.is_open()){
                    fr.close();
                    std::rename(newsjm.c_str(), orisjm.c_str());
                }
            }
        }
    }
}
