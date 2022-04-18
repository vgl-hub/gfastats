#ifndef GFASTATS_VALIDATE_H
#define GFASTATS_VALIDATE_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>
#include <vector>
#include <set>

std::string getExePath(const std::string &argv0) {
    std::string exePath = argv0.substr(0, argv0.find_last_of("/\\")+1);
    std::replace(exePath.begin(), exePath.end(), '\\', '/');
#ifdef _WIN32
    exePath += "gfastats.exe";
#else
    exePath += "gfastats";
#endif
    return exePath;
}

std::string getFileExt(const std::string& FileName) // utility to get file extension
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}

std::vector<std::string> list_dir(const char *path) {
    std::vector<std::string> list;
    struct dirent *entry;
    DIR *dir = opendir(path);

    if (dir == NULL) {
        std::cerr << "error: unable to access " << path << std::endl;
        exit(0);
    }
    while ((entry = readdir(dir)) != NULL) {
        if(entry->d_type == DT_REG) list.push_back(std::string(entry->d_name));
    }
    closedir(dir);
    return list;
}

void get_recursive(const std::string &path, std::set<std::string> &paths) {
    if(getFileExt(path) == "tst") {
        paths.insert(path);
    } else {
        std::cout << "1.1" << std::endl;
        DIR *dir = opendir(path.c_str());
        std::cout << "1.2" << std::endl;
        if(dir != NULL) {
            std::cout << "1.3" << std::endl;
            for(const auto &file : list_dir(path.c_str())) {
                get_recursive((path+"/"+file).c_str(), paths);
            }
            std::cout << "1.4" << std::endl;
        }
        std::cout << "1.5" << std::endl;
        closedir(dir);
        std::cout << "1.6" << std::endl;
    }
}

#endif // #ifndef GFASTATS_VALIDATE_H