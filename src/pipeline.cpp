#include <future>
#include <fstream>
#include "pipeline.h"

Pipeline::Pipeline(int n, int s, const std::string& fmkf, const std::string& smkf){
    goodMarkFile = smkf;
    failMarkFile = fmkf;
    pipelist.resize(n);
    for(size_t i = 0; i < pipelist.size(); ++i){
        pipelist[i].resize(s);
    }
}

Pipeline::~Pipeline(){
    for(auto& e: pipelist){
        for(auto& f: e){
            for(auto& g: f){
                if(g){
                    delete g;
                    g = NULL;
                }
            }
        }
    }
}

void Pipeline::addRunFile(RunTask* r, int s, int t){
    pipelist[s][t].push_back(r);
}

void Pipeline::prepareRerun(){
    std::string oriSJM, newSJM, bakSJM;
    for(auto& e: pipelist){
        for(auto& f: e){
            for(auto& g: f){
                remove(g->failMarkFile.c_str());
                remove(g->goodMarkFile.c_str());
                remove(g->logFile.c_str());
                auto p = g->sjmCMD.find_last_of(" ");
                oriSJM = g->sjmCMD.substr(p + 1);
                newSJM = oriSJM + ".status";
                bakSJM = newSJM + ".bak";
                remove(bakSJM.c_str());
                std::ifstream fr(newSJM);
                if(fr.is_open()){
                    fr.close();
                    if(!forceUpdateSJM){
                        std::rename(newSJM.c_str(), oriSJM.c_str());
                    }
                }
            }
        }
    }
}

int Pipeline::runPipeline(){
    std::vector<std::future<int>> fv;
    for(size_t i = 0; i < pipelist.size(); ++i){
        fv.push_back(std::async(std::launch::async, &Pipeline::runStage, this, i));
    }
    for(auto& e: fv){
        retValue += std::abs(e.get());
    }
    if(retValue){
        remove(goodMarkFile.c_str());
        std::ofstream fw(failMarkFile.c_str());
        fw.close();
    }else{
        remove(failMarkFile.c_str());
        std::ofstream fw(goodMarkFile);
        fw.close();
    }
    return retValue;
}

int Pipeline::runStage(int s){
    int ret = 0;
    std::vector<std::future<int>> fv;
    for(size_t i = 0; i < pipelist[s].size(); ++i){
        fv.clear();
        for(size_t j = 0; j < pipelist[s][i].size(); ++j){
            fv.push_back(std::async(std::launch::async, &Pipeline::runTask, this, std::ref(pipelist[s][i][j])));
        }
        for(size_t k = 0; k < fv.size(); ++k){
            pipelist[s][i][k]->retValue = fv[k].get();
            if(pipelist[s][i][k]->retValue){
                remove(pipelist[s][i][k]->goodMarkFile.c_str());
                std::ofstream fw(pipelist[s][i][k]->failMarkFile);
                fw.close();
                ret = 1;
            }else{
                remove(pipelist[s][i][k]->failMarkFile.c_str());
                std::ofstream fw(pipelist[s][i][k]->goodMarkFile);
                fw.close();
            }
        }
    }
    return ret;
}

int Pipeline::runTask(RunTask* r){
    int ret = std::system(r->sjmCMD.c_str());
    return ret;
}
