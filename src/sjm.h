#ifndef SJM_H
#define SJM_H

#include <future>
#include <string>
#include <vector>
#include <ostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <algorithm>

namespace sjm{
    // sjm specific args
    struct args{
        const int minstage = 1;
        const int maxstage = 10;
        std::string sample_list;
        std::string reg;
        std::vector<int> ana_marker;
        int ini_marker = 1;
        int end_marker = 9;
        std::string out_dir = "./";
        std::string bin_dir = "./";
        std::string db_dir = "./";
        std::string sjm_dir = "00.sjm";
        std::string cut_dir = "01.cutadapt";
        std::string spl_dir = "02.splitread";
        std::string fil_dir = "03.filtdb";
        std::string dfq_dir = "04.seqtk";
        std::string aln_dir = "05.alignment";
        std::string mkd_dir = "06.markdup";
        std::string bqc_dir = "07.bamqc";
        std::string fus_dir = "08.fusion";
        std::string exp_dir = "09.express";
        std::string rep_dir = "report";
        std::string log_dir = "log";
        std::string dfq_vol = "VOL";
        bool rerun = false;
        bool local = false;
        bool gensjm = false;
        bool showmark = false;
        std::string version = "0.0.0";
    };

    struct fmap{
        std::string fq1, fq2;
        std::string outdir, pre;
        fmap(const std::string& r1, const std::string r2, const std::string& otd, const std::string& pre) :
            fq1(r1), fq2(r2), outdir(otd), pre(pre){};
        friend std::ostream& operator<<(std::ostream& os, const fmap& f){
            os << "<Files>\n";
            os << f.fq1 << "\n" << f.fq2 << "\n\n";
            os << "<Options>\n";
            os << "MonoPath=/share/public/software/mono/mono/bin/mono\n";
            os << "PairedEnd=True\nRnaMode=True\nThreadNumber=8\nGzip=False\nOutputFusionReads=True\n\n";
            os << "<Output>\n";
            os << "TempPath=" << f.outdir << "/" << f.pre << "_tmp\n";
            os << "OutputPath=" << f.outdir << "/" << f.pre << "\n";
            os << "OutputName=" << f.pre << "\n";
            return os;
        }
    };

    // a job is a execute unit of sjm
    struct job{
        std::string begin = "job_begin";
        std::string end = "job_end";
        std::string sep = "  ";
        std::pair<std::string, std::string> name = {"name", ""};
        std::pair<std::string, std::string> memory = {"memory", ""}; // hardlimit(if use more than this amount, process will crash) in sjm, so I do not write it into sjm
        std::pair<std::string, std::string> queue = {"queue", ""};
        std::pair<std::string, std::string> slots = {"slots", ""}; // slots in PE if write it into sjm, so I use it as a convenient key to specify memory
        std::pair<std::string, std::string> cmd = {"cmd", ""};
        std::pair<std::string, std::string> envexp = {"export", ""};
        std::pair<std::string, std::string> workdir = {"directory", ""};
        std::pair<std::string, std::string> host = {"host", ""};
        std::pair<std::string, std::string> sopt = {"sched_options", "-V"};
        std::pair<std::string, std::string> status = {"status", ""};
        int stage_marker = 0;
        std::string o1, o2; // store 1/2 output file path;

        job(const int& m) : stage_marker(m){};
        job() = default;

        static std::ostream& append_item(std::ostream& os, const std::pair<std::string, std::string>& item){
            if(item.second.empty()){
                return os;
            }else{
                os << "  " << item.first << " " << item.second << "\n";
                return os;
            }
        }

        friend std::ostream& operator<<(std::ostream& os, const job& j){
            os << j.begin << "\n" 
               << "  " << j.name.first << " " << j.name.second << "\n"
               << "  " << j.cmd.first << " " << j.cmd.second << "\n";
            append_item(os, j.queue);
            append_item(os, j.envexp);
            append_item(os, j.workdir);
            append_item(os, j.status);
            append_item(os, j.sopt);
            append_item(os, j.host);
            os << j.end << "\n\n";
            return os;
        }
    };

    // a task contain a sequences of jobs
    struct task{
        std::vector<std::vector<job>> joblist;
        std::pair<std::string, std::string> logdir = {"log_dir", ""};
        friend std::ostream& operator<<(std::ostream& os, const task& t){
            for(auto& e: t.joblist){
                for(auto& f: e){
                    os << f;
                }
            }
            if(!t.logdir.second.empty()){
                os << "  " << t.logdir.first << " " << t.logdir.second << "\n";
            }
            if(t.joblist.size() <= 1){
                return os;
            }
            for(size_t i = 1; i < t.joblist.size(); ++i){
                for(auto& e: t.joblist[i-1]){
                    for(auto& f: t.joblist[i]){
                        os << "order  " << e.name.second << " before " << f.name.second << "\n";
                    }
                }
            }

            return os;
        }
    };

    // a pipeline contain a sequences of sjm files
    struct pipeline{
        std::vector<std::vector<std::vector<std::string>>> pipelist;
        std::vector<int> njob;
        int run_pipe();
        void pre_rerun();
    };

    // show stage_marker<->ana_stage relationship
    void show_mark();
    // update args after commandline args parsed
    void update_args(sjm::args& a);
    // generate subdirectories
    void gen_dir(const sjm::args& a);
    // generate fastp job
    void gen_fastp_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j);
    // generate splitr job
    void gen_splitr_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& sconf, const std::string& pre, sjm::job& j);
    // generate filtdb job
    void gen_filtdb_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j);
    // downsample fq
    void gen_seqtk_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j);
    // generate alignment job
    void gen_align_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j);
    // generate mkdup job
    void gen_mkdup_job(const sjm::args& a, const std::string& bam, const std::string& pre, sjm::job& j);
    // generate bamqc job
    void gen_bamqc_job(const sjm::args& a, const std::string& bam, const std::string& pre, sjm::job& j);
    // generate fusion job
    void gen_fusion_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j);
    // generate express job
    void gen_express_job(const sjm::args& a, const std::string& lib1, const std::string& lib2, const std::string& pre, sjm::job& j); 
    // generate report job
    void gen_report_job(const sjm::args& a, const std::string& lib);
    // generate prelib task
    void gen_prelib_task(const sjm::args& a, sjm::pipeline& p);
    // generate analib task
    void gen_analib_task(const sjm::args& a, sjm::pipeline& p);
    // ini pipeline
    void init_pipeline(const sjm::args& a, sjm::pipeline& p);
    // run sjm
    int run_sjm(const std::string& sjm);
    // run task
    int run_task(const std::vector<std::vector<std::string>>& task);
}

#endif
