/*
USAGE:
test <path to test folder or files>

EXAMPLE:
build/bin/gfastats-validate testFiles
build/bin/gfastats-validate testFiles/random1.fasta testFiles/random2.gfa2.gfa.gz


*/

#include <algorithm>
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
#include <regex>

bool verbose = false, veryVerbose = false, printCommand = false;
const std::string tmp = "tmp.txt";
const std::string err = "err.txt";
bool pass = true;

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
        fprintf(stderr, "error: unable to access <%s>\n", path);
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
        DIR *dir = opendir(path.c_str());
        if(dir != NULL) {
            for(const auto &file : list_dir(path.c_str())) {
                get_recursive((path+"/"+file).c_str(), paths);
            }
        }
        closedir(dir);
    }
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

    bool retval=true;
    bool fail;
    while(std::getline(istreamActual, line)) {
        colonIndex = line.find(':');
        if(colonIndex == std::string::npos) continue;
        std::string key = line.substr(0, colonIndex);
        std::string actual = line.substr(colonIndex+1);
        fail = expected[key] != actual;
        retval = retval && !fail;
        if(veryVerbose || fail) {
            diffs[key] = {expected[key], actual};
        }
    }

    return retval;
}

void printFAIL(const char *m1="", const char *m2="", const char *m3="", const char *m4="\n") {
    pass = false;
    printf("\033[0;31mFAIL\033[0m %s %s %s %s", m1, m2, m3, m4);
}

void printPASS(const char *m1="", const char *m2="", const char *m3="", const char *m4="\n") {
    printf("\033[0;32mPASS\033[0m %s %s %s %s", m1, m2, m3, m4);
}

int main(int argc, char **argv) {
    if (argc == 1) { // test with no arguments
        printf("gfastats-validate <path to test folder and/or files>\n");
        exit(0);
    }

    int opt;
    while((opt = getopt(argc, argv, "vVc")) != -1) 
    {
        switch(opt) 
        {
        case 'V':
            veryVerbose = true;
        case 'v':
            verbose = true;
            break;
        case 'c':
            printCommand = true;
            break;
        }
    }

    std::set<std::string> input_files;

    for(int i=1; i<argc; ++i) {
        get_recursive(argv[i], input_files);
    }

    std::string argv0 = std::string(argv[0]);
    std::string exePath = argv0.substr(0, argv0.find_last_of("/\\")+1);
    std::replace(exePath.begin(), exePath.end(), '\\', '/');

    exePath += "gfastats";

    std::string line;
    std::ifstream istream, exp, actOutput, *expOutput;
    for(const auto &input_file : input_files) {
        istream.open(input_file);
        if(!istream) {
            printFAIL(input_file.c_str(), "couldn't open test file");
            continue;
        }
        std::getline(istream, line);
        line.erase(remove(line.begin(), line.end(), '\r'), line.end());
        line.erase(remove(line.begin(), line.end(), '\n'), line.end());
        std::string cmd = exePath+line+" > "+tmp+" 2>"+err;
        if(printCommand) printf("%s\n", cmd.c_str());

        if(system(cmd.c_str()) != EXIT_SUCCESS) {
            printFAIL(input_file.c_str(), "runtime error");
            std::ifstream errfstream;
            errfstream.open(err);
            if(!errfstream) printf("    error: couldn't open err.txt\n");
            for(std::string line; std::getline(errfstream, line);) {
                printf("    %s\n", line.c_str());
            }
            errfstream.close();
            istream.close();
            continue;
        }


        std::getline(istream, line);
        exp.open(line);
        if(exp) {
            expOutput = &exp; // seperate expected output file
        } else if(line == "embedded") {
            expOutput = &istream;
        } else {
            printf("%d;;;%s;;;\n", line.length(), "");
            printFAIL("couldn't open expected output");
            continue;
        }

        actOutput.open(tmp);
        std::vector<std::pair<std::string, std::string>> diffs;
        while(!actOutput.eof() || !expOutput->eof()) {
            std::string l1, l2;
            std::getline(actOutput, l1);
            std::getline(*expOutput, l2);
            if(l1 != l2) diffs.push_back(std::pair<std::string, std::string>(l1, l2));
        }
        actOutput.close();

        exp.close();
        istream.close();

        if(diffs.size() > 0) {
            printFAIL(input_file.c_str(), "expected output did not match actual output");
            if(verbose)
            for(const auto &pair : diffs) {
                printf("    expected: %s\n      actual: %s\n", pair.second.c_str(), pair.first.c_str());
            }
            continue;
        }

        printPASS(input_file.c_str());
    } 

    if(remove(tmp.c_str()) != 0) {
        fprintf(stderr, "error deleting temp file <%s>\n", tmp.c_str());
    }

    exit(pass ? EXIT_SUCCESS : EXIT_FAILURE);
}
