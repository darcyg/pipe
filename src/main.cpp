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
    app.add_option("-s,--sample_list", a.sample_list, "sample list file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-r,--region_file", a.reg, "region bed file")->required(true)->check(CLI::ExistingFile);
    app.add_option("-v,--reads_count", a.dfq_vol, "down sample fastq reads number");
    app.add_option("-o,--out_dir", a.out_dir, "output directory");
    app.add_option("-a,--ana_mark", a.ana_marker, "analysis stage marker")->check(CLI::Range(a.minstage, a.maxstage));
    app.add_option("-i,--ini_mark", a.ini_marker, "initial analysis stage marker")->check(CLI::Range(a.minstage, a.maxstage));
    app.add_option("-e,--end_mark", a.end_marker, "end analysis stage marker")->check(CLI::Range(a.minstage, a.maxstage));
    app.add_flag("-c,--continue", a.rerun, "continue run from last failure");
    app.add_flag("-l,--loc", a.local, "Running in localhost");
    app.add_flag("-g,--gensjm", a.gensjm, "Only generate sjms, not run pipeline");
    CLI_PARSE(app, argc, argv);
    
    std::ofstream frun("./run.sh");
    for(int i = 0; i < argc; ++i){
        frun << argv[i] << " ";
    }

    sjm::update_args(a);
    sjm::pipeline p;
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
