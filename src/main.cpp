#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <climits>
#include <mutex>
#include <libgen.h>
#include "options.h"
#include "CLI.hpp"
#include "util.h"

int main(int argc, char** argv)
{
    std::string sys_cmd = std::string(argv[0]) + " -h";
    if(argc == 1){
        std::system(sys_cmd.c_str());
        sjm::show_mark();
        return 0;
    }
    std::mutex logmtx;
    util::loginfo("parsing arguments started.", logmtx);
    Options opt;
    std::string cmp_time = std::string(__TIME__) + " " + std::string(__DATE__);
    CLI::App app("program: " + std::string(basename(argv[0])) + "\nversion: " + opt.clOpt.version + "\nupdated: " + cmp_time);
    app.add_option("-s,--slist", opt.clOpt.sample_list, "sample list file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-r,--ref", opt.clOpt.ref, "reference file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-b,--bed", opt.clOpt.reg, "bed region file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-t,--gset", opt.clOpt.gset, "gene list")->check(CLI::ExistingFile);
    app.add_option("-v,--vread", opt.clOpt.dfq_vol, "fastq  subset reads number, 0 to use smallest reads num in sublibraries");
    app.add_option("-o,--out", opt.clOpt.out_dir, "output directory");
    app.add_option("-a,--amark", opt.clOpt.ana_marker, "analysis marker range")->check(CLI::Range(opt.clOpt.minstage, opt.clOpt.maxstage));
    app.add_option("-i,--imark", opt.clOpt.ini_marker, "initial analysis marker")->check(CLI::Range(opt.clOpt.minstage, opt.clOpt.maxstage));
    app.add_option("-e,--emark", opt.clOpt.end_marker, "end analysis marker")->check(CLI::Range(opt.clOpt.minstage, opt.clOpt.maxstage));
    app.add_option("-q,--queue", opt.clOpt.queue, "queue to run tasks");
    CLI::Option* prerun = app.add_flag("-c,--ctd", opt.clOpt.rerun, "continue from last failure");
    app.add_flag("-l,--loc", opt.clOpt.local, "run in localhost");
    app.add_flag("-g,--gen", opt.clOpt.gensjm, "generate sjms, not run tasks");
    app.add_flag("-u,--update", opt.clOpt.update, "update command to execute")->needs(prerun);
    CLI_PARSE(app, argc, argv);
    util::loginfo("parsing arguments finished.");
     
    util::loginfo("preparing output parent directory.");
    util::makedir(opt.ioOpt.out_dir);
    util::loginfo("output parent directory prepared.");
    
    util::loginfo("updating arguments.");
    opt.updateOptions();
    util::loginfo("arguments updated.");
    
    util::loginfo("construct pipeline object.");
    sjm::pipeline p(a.update);
    util::loginfo("pipeline object construced.");

    util::loginfo("initialize pipeline.");
    sjm::init_pipeline(a, p);
    util::loginfo("pipeline initialized.");

    util::loginfo("generate subdirectories.");
    sjm::gen_dir(a);
    util::loginfo("subdirectories generated.");

    util::loginfo("generate library prepare sub pipeline.");
    sjm::gen_prelib_task(a, p);
    util::loginfo("library prpare sub pipeline generated.");
    
    util::loginfo("generate library analysis pipeline.");
    sjm::gen_analib_task(a, p);
    util::loginfo("library analysis pipeline generated.");
    
    if(a.rerun){
        util::loginfo("resume last running of pipeline started.");
        p.pre_rerun();
        util::loginfo("finished resume running.");
    }
    if(!a.gensjm){
        util::loginfo("running pipeline now.");
        p.run_pipe();
        util::loginfo("pipeline finished.");
    }
}
