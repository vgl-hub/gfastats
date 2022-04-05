#include <gfastats-validate.h>
#include <map>

int main(int argc, char **argv) {
    std::cout << "WARNING: only run this program if gfastats is in a working state" << std::endl;
    std::cout << "continue? (Y/N) ";
    std::string input;
    std::cin >> input;
    if(input != "Y" && input != "y") {
        std::cout << "test generation cancelled" << std::endl;
        std::exit(0);
    }
    std::cout << "generating test files..." << std::endl;

    std::string exePath = getExePath(argv[0]);

    const std::map<std::set<std::string>, std::vector<std::string>> ext_args = {
        {{"fasta", "fastq"}, {"", "-b a", "-b c", "-b s", "--homopolymer-compress 5"}},
        {{"gfa"}, {""}}
    //  {{set of file extensions}, {list of command line args to run it with}}
    };

    const std::set<std::string> exclude {"agp", "sak"};

    for(const std::string &file : list_dir("testFiles")) {
        std::string ext = getFileExt(file);
        if(exclude.count(ext)) continue;

        int i = 0;
        for(auto pair : ext_args) {
            if(!pair.first.count(ext)) continue;
            for(auto args : pair.second) {
                std::string tstFile = "validateFiles/"+file+std::to_string(i)+".tst";
                std::ofstream ostream;
                ostream.open(tstFile);
                ostream << "testFiles/" << file << " " << args << "\nembedded" << std::endl;
                ostream.close();
                std::string cmd = exePath+" "+args+" testFiles/"+file+" >> "+tstFile;
                system(cmd.c_str());
                ++i;
            }
        }
    }

    std::exit(0);
}