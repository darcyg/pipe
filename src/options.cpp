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
        logfile = util::joinpath(ioOpt.fil_dir, libName + ".filter.json");
        std::ifstream fr1(logfile);
        jsn::json j;
        fr1 >> j;
        readsNum.push_back(j["FilterResult"]["Summary"]["ReadsGot"]);
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
    std::vector<std::pair<std::string, std::string>> stg;
    stg.push_back({"split read",        "spliter"});
    stg.push_back({"cutadapter and qc", "fqtool"});
    stg.push_back({"filter rrna",       "filter"});
    stg.push_back({"downsample fastq",  "seqtk"});
    stg.push_back({"genome slignment",  "bwa"});
    stg.push_back({"markdup",           "duplexer"});
    stg.push_back({"bam QC",            "bamqc"});
    stg.push_back({"fusion calling",    "fusionMap"});
    stg.push_back({"express quant",     "kallisto"});
    stg.push_back({"report",            "genrpt"});
    
    std::map<std::string, std::string> sver;
    sver["spliter"]    = "0.0.0";
    sver["fqtool"]     = "0.0.0";
    sver["filter"]     = "0.0.0";
    sver["seqtk"]      = "1.3-r106";
    sver["bwa"]        = "0.7.17-r1188";
    sver["duplexer"]   = "0.0.0";
    sver["bamqc"]      = "0.0.0";
    sver["fusionMap"]  = "10.0.1.29";
    sver["kallisto"]   = "0.45.1";
    sver["genrpt"]     = "0.0.0";
    
    std::cout << std::left;
    std::cout << "  ┌------┬-----------------┬---------┬------------┐" << std::endl;
    std::cout << "  |" << std::setw(6) << "Marker" << "|" << std::setw(17) << "Analysis" << "|";
    std::cout << std::setw(9) << "Software" << "|" << std::setw(12) << "Version" << "|" << std::endl;
    for(size_t i = 0; i < stg.size(); ++i){
    std::cout << "  |------┼-----------------┼---------┼------------|" << std::endl;
        std::cout << "  |" << std::setw(6) << i + 1 << "|" << std::setw(17) << stg[i].first << "|";
        std::cout << std::setw(9) << stg[i].second << "|" << std::setw(12) << sver[stg[i].second] << "|" << std::endl;
    }
    std::cout << "  └------┴-----------------┴---------┴------------┘" << std::endl;
}
