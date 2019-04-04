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

    std::ifstream fr(clOpt.sample_list);
    std::string line;
    int count = 0;
    while(std::getline(fr, line)){
        ++count;
    }
    nSamples = count;
    goodMarkFile = ioOpt.log_dir + "/SUCCESS";
    failMarkFile = ioOpt.log_dir + "/FAIL";
    if(clOpt.dfq_vol == "0"){
        updateMinFqVolMap();
    }
}

void Options::updateMinFqVolMap(){
    std::string line;
    std::vector<std::string> vstr;
    std::ifstream fr(clOpt.sample_list);
    std::getline(fr, line);
    while(std::getline(fr, line)){
        util::split(line, vstr, "\t");
        minFqVolMap[vstr[0]] = getMinFqVolOfOnePool(vstr[vstr.size() - 1]);
    }
}

size_t Options::getMinFqVolOfOnePool(const std::string& conf){
    if(!util::exists(conf)){
        util::error_exit("configure file \"" + conf + "\" does not existd!\n");
    }
    std::vector<size_t> readsNum;
    std::ifstream fr(conf);
    std::string logfile, line, libName;
    std::istringstream iss;
    while(std::getline(fr, line)){
        iss.clear();
        iss.str(line);
        iss >> libName;
        logfile = util::joinpath(ioOpt.fil_dir, libName + ".filt.log");
        std::ifstream fr1(logfile);
        std::string confLine;
        std::vector<std::string> vstr;
        std::getline(fr1, confLine);
        util::split(confLine, vstr, " ");
        int totReads = std::atoi(vstr[2].c_str());
        std::getline(fr1, confLine);
        util::split(confLine, vstr, " ");
        int filtReads = std::atoi(vstr[3].c_str());
        readsNum.push_back(totReads - filtReads);
        fr1.close();
    }
    size_t minReadsNum = readsNum[0];
    for(size_t i = 1; i < readsNum.size(); ++i){
        minReadsNum = std::min(minReadsNum, readsNum[i]);
    }
    fr.close();
    return minReadsNum;
}

void Options::genDirectory(){
    util::makedir(ioOpt.sjm_dir);
    util::makedir(ioOpt.cut_dir);
    util::makedir(ioOpt.spl_dir);
    util::makedir(ioOpt.fil_dir);
    util::makedir(ioOpt.dfq_dir);
    util::makedir(ioOpt.aln_dir);
    util::makedir(ioOpt.mkd_dir);
    util::makedir(ioOpt.bqc_dir);
    util::makedir(ioOpt.fus_dir);
    util::makedir(ioOpt.exp_dir);
    util::makedir(ioOpt.rep_dir);
    util::makedir(ioOpt.log_dir);
}

void Options::showMark(){
    std::vector<std::string> stg = {"cutadaptor by fastp", "split read by splitr", "filter rrna by filtdb",
                                    "downsample fastq by seqtk", "genome slignment by bwa", "markdup by mkdp",
                                    "bam QC by bamqc", "fusion calling by fusionMap", "express quant by kallisto",
                                    "report by genrep"};
    std::cout << std::left;
    std::cout << "  ┌--------------------------------------┐" << std::endl;
    std::cout << "  |" << std::setw(10) << "Marker" << std::setw(28) << "Analysis" << "|" << std::endl;
    for(size_t i = 0; i < stg.size(); ++i){
        std::cout << "  |--------------------------------------|" << std::endl;
        std::cout << "  |" << std::setw(10) << i + 1 << std::setw(28) << stg[i] << "|" << std::endl;
    }
    std::cout << "  └--------------------------------------┘" << std::endl;
}
