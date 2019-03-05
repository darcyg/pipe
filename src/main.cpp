#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "libgen.h"
#include "CLI.hpp"
#include "sjm.h"
#include "util.h"

int main(int argc, char** argv)
{
    std::string sys_cmd = std::string(argv[0]) + " -h";
    if(argc == 1){
        std::system(sys_cmd.c_str());
        return 0;
    }

    sjm::args a;
    std::string cmp_time = std::string(__TIME__) + " " + std::string(__DATE__);
    CLI::App app("program: " + std::string(basename(argv[0])) + "\nversion: " + a.version + "\nupdated: " + cmp_time);
    app.add_flag("-l,--loc", a.local, "Running in localhost");
    
    CLI::App* pgen = app.add_subcommand("gen", "gen pipeline");
    pgen->add_option("-s,--sample_list", a.sample_list, "sample list file")->required(true)->check(CLI::ExistingFile);
    pgen->add_option("-r,--region_file", a.reg, "region bed file")->required(true)->check(CLI::ExistingFile);
    pgen->add_option("-v,--reads_count", a.dfq_vol, "down sample fastq reads number");
    pgen->add_option("-o,--out_dir", a.out_dir, "output directory")->required(false);
    
    CLI::App* prun = app.add_subcommand("run", "run pipeline");
    prun->add_option("-s,--sample_list", a.sample_list, "sample list file")->required(true)->check(CLI::ExistingFile);
    prun->add_option("-r,--region_file", a.reg, "region bed file")->required(true)->check(CLI::ExistingFile);
    prun->add_option("-v,--reads_count", a.dfq_vol, "down sample fastq reads number");
    prun->add_option("-o,--out_dir", a.out_dir, "output directory");
    prun->add_flag("-c,--continue", a.rerun, "continue run from last failure");
    
    CLI_PARSE(app, argc, argv);

    a.sample_list = util::get_abspath(a.sample_list);
    a.out_dir = util::get_abspath(a.out_dir);
    a.bin_dir = util::get_dirname(util::get_abspath(argv[0]));
    a.db_dir = util::get_dirname(a.bin_dir) + "/db/";

    sjm::pipeline p;
    sjm::init_pipeline(a, p);

    if(pgen->parsed()){
        sjm::gen_dir(a);
        sjm::gen_prelib_task(a, p);
        sjm::gen_analib_task(a, p);
    }else if(prun->parsed()){
        sjm::gen_dir(a);
        sjm::gen_prelib_task(a, p);
        sjm::gen_analib_task(a, p);
        if(a.rerun){
            p.pre_rerun();
        }
        p.run_pipe();
    }else{
        return 0;
    }
}
