#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <climits>
#include "libgen.h"
#include "CLI.hpp"
#include "sjm.h"
#include "util.h"

int main(int argc, char** argv)
{
    std::string sys_cmd = std::string(argv[0]) + " -h";
    if(argc == 1){
        std::system(sys_cmd.c_str());
        sjm::show_mark();
        return 0;
    }

    sjm::args a;
    std::string cmp_time = std::string(__TIME__) + " " + std::string(__DATE__);
    CLI::App app("program: " + std::string(basename(argv[0])) + "\nversion: " + a.version + "\nupdated: " + cmp_time);
    app.add_option("-s,--slist", a.sample_list, "sample list file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-r,--ref", a.ref, "reference file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-b,--bed", a.reg, "bed region file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-v,--vread", a.dfq_vol, "fastq  subset reads number");
    app.add_option("-o,--out", a.out_dir, "output directory");
    app.add_option("-a,--amark", a.ana_marker, "analysis marker range")->check(CLI::Range(a.minstage, a.maxstage));
    app.add_option("-i,--imark", a.ini_marker, "initial analysis marker")->check(CLI::Range(a.minstage, a.maxstage));
    app.add_option("-e,--emark", a.end_marker, "end analysis marker")->check(CLI::Range(a.minstage, a.maxstage));
    app.add_option("-q,--queue", a.queue, "queue to run tasks");
    CLI::Option* prerun = app.add_flag("-c,--ctd", a.rerun, "continue from last failure");
    app.add_flag("-l,--loc", a.local, "run in localhost");
    app.add_flag("-g,--gen", a.gensjm, "generate sjms, not run tasks");
    app.add_flag("-u,--update", a.update, "update command to execute")->needs(prerun);
    CLI_PARSE(app, argc, argv);
     
    util::make_dirs(a.out_dir);
    sjm::update_args(a);
    sjm::pipeline p(a.update);
    sjm::init_pipeline(a, p);
    sjm::gen_dir(a);
    sjm::gen_prelib_task(a, p);
    sjm::gen_analib_task(a, p);
    if(a.rerun){
        p.pre_rerun();
    }
    if(!a.gensjm){
        p.run_pipe();
    }
}
