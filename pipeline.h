#ifndef PIPE_H
#define PIPE_H

#include <string>
#include <vector>

namespace pipeline{
    struct args{
        std::string sample_list;
        std::string sample_conf;
        std::string out_dir = "./";
        std::string mode;
        bool sge = true;
        std::string version = "0.0.0";
    };

    struct opts{
        std::vector<std::string> sub_dirs = {"01.cutadapt", "02.splitread", "03.alignment", "04.markdup", "05.bamqc", "06.fusion"};
        std::vector<std::vector<std::string>> sub_cmds, sub_sucfs, sub_erfs;
        int tot_cmd = 0;
        std::string bin_dir;
        std::string db_dir;
    };

    void gen_dir(const pipeline::args& args, const pipeline::opts& opts);
    void gen_sh(const pipeline::args& args, pipeline::opts& opts);
    void run_sh(const pipeline::args& args, const pipeline::opts& opts);
    int exec_cmd(const std::string& cmd, const std::string& mkf, const std::string& erf);
    void show_cmd(const pipeline::opts& opts);
}

#endif
