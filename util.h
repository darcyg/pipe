#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <string>
#include <sys/stat.h>
#include <cerrno>
#include <cstdio>

#ifdef WINDOWS
#include <direct.h>
#define get_curent_dir _getcwd
#define get_real_path _realpath
#else
#include <unistd.h>
#define get_curent_dir getcwd
#define get_real_path realpath
#endif

namespace util{
    bool is_dir_exists(const std::string& path);
    bool make_dirs(const std::string& path);
    std::string get_cwd(void);
    std::string get_abspath(const std::string& path);
    std::string get_dirname(const std::string& path);
    std::string get_basename(const std::string& path);
}

#endif
