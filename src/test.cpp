#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <vector>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <map>
#include <gfastats-functions.h>

const char *tmp = "tmp.txt";
const char *err = "err.txt";
std::string curPath, exePath, tstPath, expPath;

bool hasValidTestFileExtension(const std::string &file) {
    std::string ext = getFileExt(file);
    if(ext == "gz") {
        return hasValidTestFileExtension(rmFileExt(file));
    }
    return (ext == "fasta" || ext == "fastq" || ext == "gfa");
}

std::vector<std::string> list_dir(const char *path) {
    std::vector<std::string> list;
    struct dirent *entry;
    DIR *dir = opendir(path);

    if (dir == NULL) {
        fprintf(stderr, "error: unable to access <%s>\n", path);
        exit(0);
    }
    while ((entry = readdir(dir)) != NULL) {
        if(entry->d_type == DT_REG && hasValidTestFileExtension(entry->d_name)) list.push_back(std::string(entry->d_name));
    }
    closedir(dir);
    return list;
}

bool test(const std::string &testFile) {
    const std::string expectedOutputPath = tstPath+testFile+".tst";

    std::ifstream istreamActual, istreamExpected;
    istreamExpected.open(expectedOutputPath);
    istreamActual.open(tmp);
    if(!istreamActual || !istreamExpected) return false;

    std::map<std::string, std::string> expected;
    std::string line;
    int colonIndex;
    while(std::getline(istreamExpected, line)) {
        colonIndex = line.find(':');
        if(colonIndex == std::string::npos) continue;
        expected[line.substr(0, colonIndex)] = line.substr(colonIndex+1);
    }

    while(std::getline(istreamActual, line)) {
        colonIndex = line.find(':');
        if(colonIndex == std::string::npos) continue;
        if(expected[line.substr(0, colonIndex)] != line.substr(colonIndex+1)) {
            return false;
        }
    }

    return true;
}

int main(int argc, char **argv) {
    if (argc == 1 || argc > 4) { // test with no arguments
        printf("test <path to gfastats.exe> <path to test folder> <optional explicit extra tests file>\n");
        exit(0);
    }

    exePath = std::string(argv[1]);
    tstPath = std::string(argv[2])+'/';
    expPath = std::string(argc == 4 ? argv[3] : "");

    std::ifstream istreamTst, istreamOut;

    char cmd[100];
    for (const auto &file : list_dir(tstPath.c_str())) {
        sprintf(cmd, "%s %s%s > %s 2>%s", exePath.c_str(), tstPath.c_str(), file.c_str(), tmp, err);
        // uncomment to print the commands before running them
        // std::cout << cmd << std::endl;
        if(system(cmd) == EXIT_SUCCESS) {
            printf("%s %s\n", (test(file) ? "PASS" : "FAIL"), file.c_str());
        }
    }

    if(expPath != "") {
        std::ifstream istreamExtra;
        istreamExtra.open(expPath);
        if(istreamExtra) {
            
        }
        istreamExtra.close();
    }

    if(remove(tmp) != 0) {
        fprintf(stderr, "error deleting temp file <%s>\n", tmp);
    }

    exit(EXIT_SUCCESS);
}