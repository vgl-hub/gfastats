/*
USAGE:
test <path to test folder or files>

EXAMPLE:
build/bin/test testFiles
build/bin/test testFiles/random1.fasta testFiles/random2.gfa2.gfa.gz


*/

#include <iostream>
#include <fstream>
#include <string>
#include <dirent.h>
#include <vector>
#include <stdio.h>
#include <unistd.h>
#include <limits.h>
#include <map>
#include <set>

bool verbose = false, veryVerbose = false;
const char *tmp = "tmp.txt";
const char *err = "err.txt";
std::string curPath, exePath;

std::string rmFileExt(const std::string& path) { // utility to strip file extension from file
    if (path == "." || path == "..")
        return path;

    size_t pos = path.find_last_of("\\/.");
    if (pos != std::string::npos && path[pos] == '.')
        return path.substr(0, pos);

    return path;
}

std::string getFileExt(const std::string& FileName) // utility to get file extension
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}

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

std::map<std::string, std::pair<std::string, std::string>> diffs;

bool test(const std::string &testFile) {
    diffs.clear();
    const std::string expectedOutputPath = testFile+".tst";

    std::ifstream istreamActual, istreamExpected;
    istreamExpected.open(expectedOutputPath);
    istreamActual.open(tmp);
    if(!istreamActual || !istreamExpected) return false;

    std::map<std::string, std::string> expected;
    std::string line;
    size_t colonIndex;
    while(std::getline(istreamExpected, line)) {
        colonIndex = line.find(':');
        if(colonIndex == std::string::npos) continue;
        expected[line.substr(0, colonIndex)] = line.substr(colonIndex+1);
    }

    while(std::getline(istreamActual, line)) {
        colonIndex = line.find(':');
        if(colonIndex == std::string::npos) continue;
        std::string key = line.substr(0, colonIndex);
        std::string actual = line.substr(colonIndex+1);
        if(expected[key] != actual) {
            diffs[key] = {expected[key], actual};
        }
    }

    return diffs.size() == 0;
}

void printFile(const char *fname) {
    std::ifstream istream;
    istream.open(fname);
    if(!istream) return;
    for(std::string line; std::getline(istream, line);) {
        printf("    %s\n", line.c_str());
    }
    istream.close();
}

int main(int argc, char **argv) {
    if (argc == 1) { // test with no arguments
        printf("test <path to test folder and/or files>\n");
        exit(0);
    }

    int opt;
    while((opt = getopt(argc, argv, "vV")) != -1) 
    {
        switch(opt) 
        {
        case 'V':
            veryVerbose = true;
        case 'v':
            verbose = true;
            break;
        }
    }

    std::map<std::string, bool> input;
    std::set<std::string> tested;

    for(int i=1; i<argc; ++i) {
        DIR *dir = opendir(argv[i]);
        if(dir != NULL || hasValidTestFileExtension(std::string(argv[i]))) input[argv[i]] = (dir != NULL);
        closedir(dir);
    }

    std::string argv0 = std::string(argv[0]);
    exePath = argv0.substr(0, argv0.find_last_of('\\')+1)+"gfastats.exe";

    char cmd[100];

    auto lambda = [&cmd, &tested](std::string file) {
        if(tested.count(file)) return;
        tested.insert(file);
        sprintf(cmd, "%s %s > %s 2>%s", exePath.c_str(), file.c_str(), tmp, err);
        if(veryVerbose) printf("%s\n", cmd);
        bool exitSuccess = system(cmd) == EXIT_SUCCESS;
        bool pass = (exitSuccess && test(file));
        printf("%s %s\n", pass ? "PASS" : "FAIL", file.c_str());
        if(veryVerbose || (!pass && verbose)) {
            if(!exitSuccess) {
                printFile(err);
            }
            else if(diffs.size() != 0) {
                for(auto const &diff : diffs)
                    printf("    key: %s; expected: %s; actual: %s\n", diff.first.c_str(), diff.second.first.c_str(), diff.second.second.c_str());
            }
            else if(!pass) {
                printf("    failed to open expected output file and/or tmp file");
            }
        }
    };

    for(const auto &i : input) {
        if(i.second) {
            for (const auto &file : list_dir(i.first.c_str())) {
                lambda(i.first+"/"+file);
            }
        } else {
            lambda(i.first);
        }
    }

    if(remove(tmp) != 0) {
        fprintf(stderr, "error deleting temp file <%s>\n", tmp);
    }

    exit(EXIT_SUCCESS);
}