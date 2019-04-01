#include "job.h"

Job::Job(){
}

Job::~Job(){
}

Job::Job(const std::string& jname, const std::string& jout, const std::string& jpre, const int& jstage){
    stage_marker = jstage;
    name.second = jname;
    workdir.second = jout + "/";
    std::string sgee = workdir.second + jpre + ".sub.e";
    std::string sgeo = workdir.second + jpre + ".sub.o";
    sopt.second.append(" -e " + sgee + " -o " + sgeo);
}

std::ostream& Job::appendItem(std::ostream& os, const std::pair<std::string, std::string>& item){
    if(item.second.empty()){
        return os;
    }else{
        os << "  " << item.first << " " << item.second << "\n";
        return os;
    }
}

std::ostream& operator<<(std::ostream& os, const Job& j){
    os << j.begin << "\n";
    os << "  " << j.name.first << " " << j.name.second << "\n";
    os << "  " << j.cmd.first << " " << j.cmd.second << "\n";
    j.appendItem(os, j.queue);
    j.appendItem(os, j.envexp);
    j.appendItem(os, j.workdir);
    j.appendItem(os, j.status);
    j.appendItem(os, j.sopt);
    j.appendItem(os, j.host);
    os << j.end << "\n\n";
    return os;
}