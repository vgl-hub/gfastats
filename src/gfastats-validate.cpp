/*
USAGE:
test <path to test folder or files>

EXAMPLE:
build/bin/gfastats-validate validateFiles
build/bin/gfastats-validate validateFiles/random1.fasta0.tst


*/

#include <gfastats-validate.h>
#include <algorithm>
#include <fstream>
#include <string>
#include <dirent.h>
#include <vector>
#include <unistd.h>
#include <limits.h>
#include <map>
#include <set>
#include <regex>

bool printCommand = false;
const std::string tmp = "tmp.txt";
const std::string err = "err.txt";
bool pass = true;

void printFAIL(const char *m1="", const char *m2="", const char *m3="", const char *m4="") {
    pass = false;
    std::cout << "\033[0;31mFAIL\033[0m " << m1 << " " << m2 << " " << m3 << " " << m4 << std::endl;
}

void printPASS(const char *m1="", const char *m2="", const char *m3="", const char *m4="") {
    std::cout << "\033[0;32mPASS\033[0m " << m1 << " " << m2 << " " << m3 << " " << m4 << std::endl;
}

bool getline(std::ifstream &s, std::string &line) {
    if(!std::getline(s, line)) return false;
    line.erase(remove(line.begin(), line.end(), '\r'), line.end());
    line.erase(remove(line.begin(), line.end(), '\n'), line.end());
    return true;
}

int main(int argc, char **argv) {
    if (argc == 1) { // test with no arguments
        std::cout << "gfastats-validate <path to test folder and/or files>" << std::endl;
        exit(EXIT_SUCCESS);
    }

    int opt;
    while((opt = getopt(argc, argv, "c")) != -1) 
    {
        switch(opt) 
        {
        case 'c':
            printCommand = true;
            break;
        }
    }

    std::set<std::string> input_files;

    for(int i=1; i<argc; ++i) {
        get_recursive(argv[i], input_files);
    }

    std::string exePath = getExePath(argv[0]);

    std::string line;
    std::ifstream istream, exp, actOutput, *expOutput;
    for(const auto &input_file : input_files) {
        istream.open(input_file);
        if(!istream) {
            printFAIL(input_file.c_str(), "couldn't open test file");
            continue;
        }
        getline(istream, line);
#ifdef _WIN32
        std::string cmd = "\"\""+exePath+"\""+" "+line+" > "+tmp+" 2>"+err+"\"";
#else
        std::string cmd = "\""+exePath+"\""+" "+line+" > "+tmp+" 2>"+err;
#endif
        if(printCommand) std::cout << cmd << std::endl;

        if(system(cmd.c_str()) != EXIT_SUCCESS) {
            printFAIL(input_file.c_str(), "runtime error");
            istream.close();
            std::ifstream errfstream;
            errfstream.open(err);
            if(!errfstream) {
                std::cout << "    error: couldn't open err.txt" << std::endl;
                continue;
            }
            for(std::string line; getline(errfstream, line);) {
                std::cout << "    " << line.c_str() << std::endl;
            }
            errfstream.close();
            continue;
        }


        getline(istream, line);
        exp.open(line);
        if(exp) {
            expOutput = &exp; // seperate expected output file
        } else if(line == "embedded") {
            expOutput = &istream;
        } else {
            printFAIL("couldn't open expected output");
            continue;
        }

        actOutput.open(tmp);
        std::string line;
        getline(*expOutput, line);
        if(line == "+++Summary+++: ") {
            getline(actOutput, line);
            std::set<std::string> exp_summary, act_summary;
            while(!actOutput.eof()) {
                getline(actOutput, line);
                act_summary.insert(line);
            }
            while(!expOutput->eof()) {
                getline(*expOutput, line);
                exp_summary.insert(line);
            }
            std::set<std::string> additions, missings;
            for(const auto &entry : exp_summary) {
                if(act_summary.count(entry) == 0) {
                    missings.insert(entry);
                }
            }
            for(const auto &entry : act_summary) {
                if(exp_summary.count(entry) == 0) {
                    additions.insert(entry);
                }
            }

            actOutput.close();
            exp.close();
            istream.close();

            if(additions.size() > 0 || missings.size() > 0) {
                printFAIL(input_file.c_str(), "expected output did not match actual output");
                std::cout << "additions:" << std::endl;
                for(const auto &addition : additions) {
                    std::cout << addition << std::endl;
                }
                std::cout << "missing:" << std::endl;
                for(const auto &missing : missings) {
                    std::cout << missing << std::endl;
                }

                continue; // to next validate file
            }
        }
        else {
            std::vector<std::pair<std::string, std::string>> diffs;

            std::string l1, l2;
            getline(actOutput, l1);
            l2 = line;
            if(l1 != l2) diffs.push_back(std::pair<std::string, std::string>(l1, l2));

            while(!actOutput.eof() || !expOutput->eof()) {
                getline(actOutput, l1);
                getline(*expOutput, l2);
                if(l1 != l2) diffs.push_back(std::pair<std::string, std::string>(l1, l2));
            }

            actOutput.close();
            exp.close();
            istream.close();

            if(diffs.size() > 0) {
                printFAIL(input_file.c_str(), "expected output did not match actual output");
                for(const auto &pair : diffs) {
                    std::cout << "    expected: " << pair.second.c_str() << std::endl << "      actual: " << pair.first.c_str() << std::endl;
                }
                continue;
            }
        }

        printPASS(input_file.c_str());
    }

    if(input_files.size() != 0 && remove(tmp.c_str()) != 0) {
        std::cerr << "error deleting temp file " << tmp.c_str() << std::endl;
    }

    exit(pass ? EXIT_SUCCESS : EXIT_FAILURE);
}
