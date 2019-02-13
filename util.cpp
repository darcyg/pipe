#include "util.h"

bool util::is_dir_exists(const std::string& path){
#if defined(_WIN32)
    struct _stat info;
    if(_stat(path.c_str(), &info) != 0){
        return false;
    }
    return (info.st_mode & _S_IFDIR) != 0;
#else
    struct stat info;
    if(stat(path.c_str(), &info) != 0){
        return false;
    }
    return (info.st_mode & S_IFDIR) != 0;
#endif
}

bool util::make_dirs(const std::string& path){
    std::string clean_path = "";
    char las = '\0';
    for(auto& e: path){
        if(e == '/' && las == '/'){
            continue;
        }else{
            clean_path += e;
            las = e;
        }
    }
#if defined(_WIN32)
    int ret = _mkdir(clean_path.c_str());
#else
    mode_t mode = 0755;
    int ret = mkdir(clean_path.c_str(), mode);
#endif
    if(ret == 0){
        return true;
    }

    switch(errno){
        case ENOENT:
            {
                int pos = clean_path.find_last_of('/');
                if(pos == std::string::npos)
#if defined(_WIN32)
                pos = clean_path.find_last_of('\\');
                if(pos == std::string::npos)
#endif
                    return false;
                if(!util::make_dirs(clean_path.substr(0, pos))){
                    return false;
                }
            }
#if defined(_WIN32)
            return 0 == _mkdir(clean_path.c_str());
#else
            return 0 == mkdir(clean_path.c_str(), mode);
#endif
        case EEXIST:
            return is_dir_exists(clean_path);
        default:
            return false;
    }
}

std::string util::get_cwd(void){
    char cpath[FILENAME_MAX];
    get_curent_dir(cpath, FILENAME_MAX);
    return cpath;
}

std::string util::get_dirname(const std::string& path){
    char cpath[FILENAME_MAX];
    get_real_path(path.c_str(), cpath);
    std::string clean_path = std::string(cpath);
    size_t found = clean_path.find_last_of("/\\");
    if(found != std::string::npos){
        return clean_path.substr(0, found);
    }else{
        return clean_path;
    }
}

std::string util::get_basename(const std::string& path){
    char cpath[FILENAME_MAX];
    get_real_path(path.c_str(), cpath);
    std::string clean_path = std::string(cpath);
    size_t found = clean_path.find_last_of("/\\");
    if(found != std::string::npos){
        return clean_path.substr(found + 1, clean_path.size());
    }else{
        return clean_path;
    }
}

std::string util::get_abspath(const std::string& path){
    char cpath[FILENAME_MAX];
    get_real_path(path.c_str(), cpath);
    std::string abspath = std::string(cpath);
    return abspath;
}
