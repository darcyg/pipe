#include "pipeline.h"
#include "util.h"
#include "ThreadPool.h"
#include <fstream>
#include <sstream>

void pipeline::gen_dir(const pipeline::args& args, const pipeline::opts& opts){
    std::string sep = "/";
#if defined(_WIN32)
    sep = "\\";
#endif
    std::string temp_path;
    for(auto& e: opts.sub_dirs){
        temp_path = args.out_dir + sep + e;
        util::make_dirs(temp_path);
    }
}

void pipeline::gen_sh(const pipeline::args& args, pipeline::opts& opts){
    opts.sub_cmds.resize(opts.sub_dirs.size());
    opts.sub_sucfs.resize(opts.sub_dirs.size());
    opts.sub_erfs.resize(opts.sub_dirs.size());
    std::string samplename, flowcell, libname, read1, read2, tmpstr;
    std::istringstream iss;
    std::ofstream fw;
    std::ifstream fr(args.sample_list);
    std::getline(fr, tmpstr);
    while(std::getline(fr, tmpstr)){
        iss.str(tmpstr);
        iss >> samplename >> flowcell >> libname >> read1 >> read2;
        //generate cutadaptor scripts
        std::string cut_dir = args.out_dir + "/" + opts.sub_dirs[0] + "/";
        std::string suc_f = cut_dir + libname + ".cut.SUCCESS";
        std::string cmd_str = opts.bin_dir + "/fastp --detect_adapter_for_pe";
        cmd_str += " -i " + read1 + " -I " + read2;
        cmd_str += " -o " + cut_dir + libname + ".R1.fq.gz -O " + cut_dir + libname + ".R2.fq.gz >";
        cmd_str += cut_dir + libname + ".cut.log 2>&1 && touch ";
        cmd_str += suc_f + "\n";
        std::string sh_file = cut_dir + libname + ".cut.sh";
        fw.open(sh_file.c_str());
        fw << cmd_str;
        fw.close();
        if(args.sge){
            cmd_str = "qsub -e " + cut_dir + libname + ".cut.e" + " -o " + cut_dir + libname + ".cut.o";
            cmd_str += " -l vf=1g -l p=1 " + sh_file;
        }
        opts.sub_cmds[0].push_back(cmd_str);
        opts.sub_sucfs[0].push_back(suc_f);
        opts.sub_erfs[0].push_back(cut_dir + libname + ".cut.e");
        ++opts.tot_cmd;
        //generate splitread scripts
        std::string spl_dir = args.out_dir + "/" + opts.sub_dirs[1] + "/";
        suc_f = spl_dir + libname + ".spl.SUCCESS";
        cmd_str = opts.bin_dir + "/splitr -b " + opts.db_dir + "/barcode/barcode.conf ";
        cmd_str += "-s " + args.sample_conf + " -r " + cut_dir + libname + ".R1.fq.gz -R ";
        cmd_str += cut_dir + libname + ".R2.fq.gz -o " + spl_dir + libname + "/ >";
        cmd_str += spl_dir + libname + ".spl.log 2>&1 && touch ";
        cmd_str += suc_f + "\n";
        sh_file = spl_dir + libname + ".spl.sh";
        fw.open(sh_file);
        fw << cmd_str;
        fw.close();
        if(args.sge){
            cmd_str = "qsub -e " + spl_dir + libname + ".spl.e" + " -o " + spl_dir + libname + ".spl.o";
            cmd_str += " -l vf=1g -l p=1 " + sh_file;
        }
        opts.sub_cmds[1].push_back(cmd_str);
        opts.sub_sucfs[1].push_back(suc_f);
        opts.sub_erfs[1].push_back(spl_dir + libname + ".spl.e");
        ++opts.tot_cmd;
        //generate alignment/bamqc/fusion scripts
        std::ifstream cfr(args.sample_conf);
        std::istringstream cis;
        std::string fina_lib, fina_tmp;
        while(std::getline(cfr, fina_tmp)){
            cis.str(fina_tmp);
            cis >> fina_lib;
            //alignment
            std::string aln_dir = args.out_dir + "/" + opts.sub_dirs[2] + "/";
            suc_f = aln_dir + fina_lib + ".aln.sort.SUCCESS";
            cmd_str = opts.bin_dir + "/bwa mem -C -t 8 " + opts.db_dir + "/refMrna/refMrna.fa ";
            cmd_str += spl_dir + libname + "/" + fina_lib + ".R1.fq " + spl_dir + libname + "/" + fina_lib + ".R2.fq 2> " + aln_dir + fina_lib + ".aln.log | ";
            cmd_str += opts.bin_dir + "/samtools view -@ 8 -o " + aln_dir + fina_lib  + ".bam >>";
            cmd_str += aln_dir + fina_lib + ".aln.log 2>&1 && ";
            cmd_str += opts.bin_dir + "/samtools sort -@ 8 -o " + aln_dir + fina_lib  + ".sort.bam " + aln_dir + fina_lib  + ".bam >";
            cmd_str += aln_dir + fina_lib + ".sort.log 2>&1 && touch ";
            cmd_str += suc_f + "\n";
            sh_file = aln_dir + fina_lib + ".aln.sh";
            fw.open(sh_file);
            fw << cmd_str;
            fw.close();
            if(args.sge){
                cmd_str = "qsub -e " + aln_dir + fina_lib + ".aln.e" + " -o " + aln_dir + fina_lib + ".aln.o";
                cmd_str += " -l vf=8g -l p=8 " + sh_file;
            }
            opts.sub_cmds[2].push_back(cmd_str);
            opts.sub_sucfs[2].push_back(suc_f);
            opts.sub_erfs[2].push_back(aln_dir + fina_lib + ".aln.e");
            ++opts.tot_cmd;
            //markdup
            std::string mkd_dir = args.out_dir + "/" + opts.sub_dirs[3] + "/";
            suc_f = mkd_dir + fina_lib + ".mkdp.SUCCESS";
            cmd_str = opts.bin_dir + "/mkdup -i " + aln_dir + fina_lib  + ".sort.bam -o ";
            cmd_str += mkd_dir + fina_lib  + ".sort.mkdup.bam >" + mkd_dir + fina_lib + ".sort.mkdup.log 2>&1 && touch ";
            cmd_str += suc_f + "\n";
            sh_file = mkd_dir + fina_lib + ".mkdp.sh";
            fw.open(sh_file);
            fw << cmd_str;
            fw.close();
            if(args.sge){
                cmd_str = "qsub -e " + mkd_dir + fina_lib + ".mkdp.e" + " -o " + mkd_dir + fina_lib + ".mkdp.o";
                cmd_str += " -l vf=1g -l p=1 " + sh_file;
            }
            opts.sub_cmds[3].push_back(cmd_str);
            opts.sub_sucfs[3].push_back(suc_f);
            opts.sub_sucfs[3].push_back(mkd_dir + fina_lib + ".mkdp.e");
            ++opts.tot_cmd;
            //bamqc
            std::string bqc_dir = args.out_dir + "/" + opts.sub_dirs[4] + "/";
            suc_f = bqc_dir + fina_lib + ".qc.SUCCESS";
            cmd_str = opts.bin_dir + "/bamqc -i " + mkd_dir + fina_lib  + ".sort.mkdup.bam -b ";
            cmd_str += opts.db_dir + "/region/reg.bed -o " + bqc_dir + " -p " + fina_lib + " -n " + fina_lib;
            cmd_str += " >" + bqc_dir + fina_lib + ".qc.log 2>&1 && touch ";
            cmd_str += suc_f + "\n";
            sh_file = bqc_dir + fina_lib + ".qc.sh";
            fw.open(sh_file);
            fw << cmd_str;
            fw.close();
            if(args.sge){
                cmd_str = "qsub -e " + bqc_dir + fina_lib + ".qc.e" + " -o " + bqc_dir + fina_lib + ".qc.o";
                cmd_str += " -l vf=1g -l p=1 " + sh_file;
            }
            opts.sub_cmds[4].push_back(cmd_str);
            opts.sub_sucfs[4].push_back(suc_f);
            opts.sub_erfs[4].push_back(bqc_dir + fina_lib + ".qc.e");
            ++opts.tot_cmd;
        }
    }
}

void pipeline::run_sh(const pipeline::args& args, const pipeline::opts& opts){
    std::vector<std::future<int>> fv;
    int sum = 0;
    for(size_t i = 0; i < opts.sub_cmds.size(); ++i){
        sum = 0;
        fv.clear();
        ThreadPool::ThreadPool pool(opts.sub_cmds[i].size());
        for(size_t j = 0; j < opts.sub_cmds[i].size(); ++j){
            fv.push_back(pool.enqueue(pipeline::exec_cmd, std::ref(opts.sub_cmds[i][j]), 
                        std::ref(opts.sub_sucfs[i][j]), std::ref(opts.sub_erfs[i][j])));
        }
        for(auto& e : fv){
            sum += e.get();
        }
        if(sum){
            break;
        }
    }
}

int pipeline::exec_cmd(const std::string& cmd, const std::string& mkf, const std::string& erf){
    int ret = std::system(cmd.c_str());
    std::ifstream fr(erf);
    std::ifstream fs(mkf);
    if(fs.is_open()){
        return 0;
    }
    while(1){
        std::this_thread::sleep_for(std::chrono::seconds(5));
        if(ret){
            return ret;
        }
        if(fr.is_open() && fr.peek() != std::ifstream::traits_type::eof()){
            return 1;
        }
        if(fs.is_open()){
            return 0;
        }
    }
}

void pipeline::show_cmd(const pipeline::opts& opts){
    for(auto& e : opts.sub_cmds){
        for(auto& f : e){
            std::cout << f << std::endl;
        }
    }
}
