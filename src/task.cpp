#include "task.h"

Task::Task(int j){
    joblist.resize(j);
}

Task::~Task(){
    for(int i = 0; i < joblist.size(); ++i){
        for(int j = 0; j < joblist[i].size(); ++j){
            if(joblist[i][j]){
                delete joblist[i][j];
                joblist[i][j] = NULL;
            }
        }
    }
}

void Task::addJob(Job* j, int i){
    joblist[i].push_back(j);
}

std::ostream& operator<<(std::ostream& os, const Task& t){
    for(auto& e: t.joblist){
        for(auto& f: e){
            os << (*f);
        }
    }
    if(!t.logdir.second.empty()){
        os << "  " << t.logdir.first << " " << t.logdir.second << "\n";
    }
    if(t.joblist.size() <= 1){
        return os;
    }
    for(size_t i = 1; i < t.joblist.size(); ++i){
        for(auto& e: t.joblist[i-1]){
            for(auto& f: t.joblist[i]){
                os << "order  " << e->name.second << " before " << f->name.second << "\n";
            }
        }
    }
    return os;
 }

