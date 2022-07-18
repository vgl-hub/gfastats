#ifndef GFASTATS_VALIDATE_H
#define GFASTATS_VALIDATE_H

#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>
#include <vector>
#include <set>
#include <string>

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

std::string rmFileExt(const std::string path) { // utility to strip file extension from file
    if (path == "." || path == "..")
        return path;

    size_t pos = path.find_last_of("\\/.");
    if (pos != std::string::npos && path[pos] == '.')
        return path.substr(0, pos);

    return path;
}

std::string getFileExt(std::string fileName) // utility to get file extension
{
    if(fileName.find_last_of(".") != std::string::npos) {
        
        if(fileName.substr(fileName.find_last_of(".")+1) == "gz") {
            
            fileName = rmFileExt(fileName);
            
            return getFileExt(fileName) + ".gz";
            
        }
        
        return fileName.substr(fileName.find_last_of(".")+1);
    }
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
        DIR *f = opendir((std::string(path)+"/"+entry->d_name).c_str());
        if(f == NULL) /*not a directory*/ list.push_back(std::string(entry->d_name));
        else closedir(f);
    }
    closedir(dir);
    return list;
}

void get_recursive(const std::string &path, std::set<std::string> &paths) {
    if(getFileExt(path) == "tst") {
        paths.insert(path);
    } else {
        DIR *dir = opendir(path.c_str());
        if(dir != NULL) {
            for(const auto &file : list_dir(path.c_str())) {
                get_recursive((path+"/"+file).c_str(), paths);
            }
            closedir(dir);
        }
    }
}

int i = 0;

void genTest(std::string exePath, const std::string &file, const std::string &args){
    std::string tstFile = "validateFiles/"+file+"."+std::to_string(i)+".tst";
    std::cout << "generating: " << tstFile << std::endl;
    std::ofstream ostream;
    ostream.open(tstFile);
    ostream << "testFiles/" << file << " " << args << "\nembedded" << std::endl;
    ostream.close();
#ifdef _WIN32
    std::string cmd = "\"\""+exePath+"\" "+args+" testFiles/"+file+" >> "+tstFile+"\"";
#else
    std::string cmd = "\""+exePath+"\" "+args+" testFiles/"+file+" >> "+tstFile;
#endif
    int exit = system(cmd.c_str());
    if (exit == EXIT_SUCCESS) {
        ostream << cmd << std::endl;
        ostream << "Command executed.";
    }
    ++i;
};

#endif // #ifndef GFASTATS_VALIDATE_H
