#include <stdlib.h>
#include <map>
#include <cstdio>

#include "validate.h"

int main(int argc, char **argv) {
    std::cout << "WARNING: only run this program if gfastats is in a working state" << std::endl;
    std::cout << "WARNING: previous validate files will be deleted" << std::endl;
    std::cout << "continue? (Y/N) ";
    std::string input;
    std::cin >> input;
    if(input != "Y" && input != "y") {
        std::cout << "validate generation cancelled" << std::endl;
        std::exit(0);
    }
    std::cout << "deleting old validate files..." << std::endl;

    for(auto &file : list_dir("validateFiles")) {
        if(getFileExt(file) != "tst") continue; // dont delete README
        file = "validateFiles/"+file;
        if(remove(file.c_str()) != 0) {
            std::cerr << "error deleting <" << file << ">" << std::endl;
            return -1;
        }
    }

    std::cout << "generating new validate files..." << std::endl;

    std::string exePath = getExePath(argv[0]);

    const std::map<std::set<std::string>, std::vector<std::string>> ext_args = {
        {{"fasta", "fasta.gz", "fastq", "fastq.gz"}, {"", "-b a", "-b c", "-b s", "--homopolymer-compress 1 -ofa"}},
        {{"gfa", "gfa.gz", "gfa2", "gfa2.gz"}, {"-o gfa2", "-o gfa", "-o fasta"}}
    //  {{set of test file extensions}, {list of command line args to run with}}
    };

    const std::map<std::set<std::string>, std::vector<std::string>> file_args = {
        {{"random1.fasta", "random1.fasta.gz", "random1.fastq", "random1.fastq.gz", "random1.gfa"}, {"-a testFiles/random1.agp --stats", "-a testFiles/random1.agp --stats -ofa"}},
        {{"random1.fasta"}, {"-k testFiles/random1.instructions.sak", "-ofa -k testFiles/random1.instructions.sak", "-ofa -k testFiles/random1.hc.sak", "-ofa -k testFiles/random1.hdc.sak"}},
        {{"random2.noseq.gfa"}, {""}},
        {{"random1.gfa2"}, {"-k testFiles/random1.gfa2.instructions.sak"}},
        {{"random1.fasta"}, {"Header2", "-r testFiles/random1.fastq"}},
        {{"random4.fasta"}, {""}},
        {{"random1.fasta.gz"}, {"-ofa"}}
        
    //  {{set of test file paths}, {list of command line args to run with}}
    };

    const std::set<std::string> exclude {"agp", "sak"};

    for(const std::string &file : list_dir("testFiles")) {
        std::string ext = getFileExt(file);
        if(exclude.count(ext)) continue;
        for(auto pair : ext_args) {
            if(!pair.first.count(ext)) continue;
            for(auto args : pair.second) {
                genTest(exePath, file, args);
            }
        }
    }

    std::fstream fstream;
    for(const auto &pair : file_args) {
        for(const std::string &file : pair.first) {
            fstream.open("testFiles/"+file);
            if(!fstream) continue;
            fstream.close();
            for(const std::string &args : pair.second) {
                genTest(exePath, file, args);
            }
        }
    }

    std::exit(EXIT_SUCCESS);
}
