#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include "libgen.h"
#include "CLI.hpp"
#include "pipeline.h"
#include "util.h"

int main(int argc, char** argv)
{
    std::string sys_cmd = std::string(argv[0]) + " -h";
    if(argc == 1){
        std::system(sys_cmd.c_str());
        return 0;
    }

    pipeline::args args;
    std::string cmp_time = std::string(__TIME__) + " " + std::string(__DATE__);
    CLI::App app("program: " + std::string(basename(argv[0])) + "\nversion: " + args.version + "\nupdated: " + cmp_time);
    app.add_flag("-s,--sge", args.sge, "Not Running in SGE");
    
    CLI::App* pgen = app.add_subcommand("gen", "generate directories/scripts");
    pgen->add_option("-l,--sample_list", args.sample_list, "sample list file")->required(true)->check(CLI::ExistingFile);
    pgen->add_option("-c,--sample_conf", args.sample_conf, "sample conf file")->required(true)->check(CLI::ExistingFile);
    pgen->add_option("-o,--out_dir", args.out_dir, "output directory")->required(false);
    
    CLI::App* prun = app.add_subcommand("run", "run pipelineline");
    prun->add_option("-l,--sample_list", args.sample_list, "sample list file")->required(true)->check(CLI::ExistingFile);
    prun->add_option("-c,--sample_conf", args.sample_conf, "sample conf file")->required(true)->check(CLI::ExistingFile);
    prun->add_option("-o,--out_dir", args.out_dir, "output directory")->required(false);
    
    CLI_PARSE(app, argc, argv);

    args.sample_list = util::get_abspath(args.sample_list);
    args.sample_conf = util::get_abspath(args.sample_conf);
    args.out_dir = util::get_abspath(args.out_dir);

    pipeline::opts opts;
    opts.bin_dir = util::get_dirname(argv[0]);
    opts.db_dir = util::get_dirname(opts.bin_dir) + "/db/";

    if(pgen->parsed()){
        pipeline::gen_dir(args, opts);
        pipeline::gen_sh(args, opts);
    }else if(prun->parsed()){
        pipeline::gen_dir(args, opts);
        pipeline::gen_sh(args, opts);
        pipeline::run_sh(args, opts);
    }else{
        return 0;
    }
}
