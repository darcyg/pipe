#include "util.h"
#include "dirutil.h"
#include "options.h"

Options::Options(){
}

Options::~Options(){
}

void Options::updateOptions(){
    clOpt.sample_list = util::abspath(clOpt.sample_list);
    if(clOpt.ana_marker.empty()){
        for(int i = clOpt.ini_marker; i <= clOpt.end_marker; ++i){
            clOpt.ana_marker.push_back(i);
        }
    }else{
        std::remove_if(clOpt.ana_marker.begin(), clOpt.ana_marker.end(), 
                       [&](int& e){return e < clOpt.ini_marker || e > clOpt.end_marker;}
                       );
    }
    if(clOpt.gset.empty()){
        clOpt.gset = ioOpt.db_dir + "gset/kgset.tsv";
    }

    ioOpt.bin_dir = util::dirname(dirutil::getExecutablePath());
    ioOpt.out_dir = util::abspath(ioOpt.out_dir);
    ioOpt.db_dir = util::dirname(ioOpt.bin_dir) + "/db/";
    
    std::string sep = "/";
    ioOpt.sjm_dir = ioOpt.out_dir + sep + ioOpt.sjm_dir;
    ioOpt.cut_dir = ioOpt.out_dir + sep + ioOpt.cut_dir;
    ioOpt.spl_dir = ioOpt.out_dir + sep + ioOpt.spl_dir;
    ioOpt.fil_dir = ioOpt.out_dir + sep + ioOpt.fil_dir;
    ioOpt.dfq_dir = ioOpt.out_dir + sep + ioOpt.dfq_dir;
    ioOpt.aln_dir = ioOpt.out_dir + sep + ioOpt.aln_dir;
    ioOpt.mkd_dir = ioOpt.out_dir + sep + ioOpt.mkd_dir;
    ioOpt.bqc_dir = ioOpt.out_dir + sep + ioOpt.bqc_dir;
    ioOpt.fus_dir = ioOpt.out_dir + sep + ioOpt.fus_dir;
    ioOpt.exp_dir = ioOpt.out_dir + sep + ioOpt.exp_dir;
    ioOpt.rep_dir = ioOpt.out_dir + sep + ioOpt.rep_dir;
    ioOpt.log_dir = ioOpt.out_dir + sep + ioOpt.log_dir;
}
