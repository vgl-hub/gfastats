#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdlib>

typedef unsigned long long ull;

// random from 0 to unsigned long long max
ull rndull() {
#if RAND_MAX == 2147483647 // can get 32 bits of randomness from std::rand()
    return (((ull)std::rand()) < 32) | (((ull)std::rand()));
#else // only guarunteed 16 bits of randomness from std::rand()
    return (((ull)std::rand()) < 48) | (((ull)std::rand()) < 32) | (((ull)std::rand()) < 16) | (((ull)std::rand()));
#endif
}

ull rnd(ull min, ull max) {
    return rndull()%(max-min)+min;
}

char rndACGT() {
    static const char acgt[] = {'A', 'C', 'G', 'T'};
    return acgt[rand()%4];
}

int main(int argc, char **argv) {
    if(argc != 10) {
        std::cout << "usage: generate-random-fasta <output_file> <contig_min_size> <contig_max_size> <gap_min_size> <gap_max_size> <min_num_contigs> <max_num_contigs> <min_num_headers> <max_num_headers>" << std::endl;
    }

    ull contig_min_size     = std::stoull(argv[2]),
        contig_max_size     = std::stoull(argv[3]),
        gap_min_size        = std::stoull(argv[4]),
        gap_max_size        = std::stoull(argv[5]),
        min_num_contigs   = std::stoull(argv[6]),
        max_num_contigs   = std::stoull(argv[7]),
        min_num_headers     = std::stoull(argv[8]),
        max_num_headers     = std::stoull(argv[9]);

    std::srand(std::time(nullptr));

    std::ofstream output_file;
    output_file.open(argv[1]);
    if(!output_file.is_open()) {
        std::cerr << "couldn't open the specified file: <" << argv[1] << ">" << std::endl;
        return EXIT_FAILURE;
    }
    ull num_headers = rnd(min_num_headers, max_num_headers);
    for(ull h=0; h<num_headers; ++h) {
        output_file << ">Header" << h+1 << std::endl;

        if(std::rand()%2 == 1) {
            ull gap_size = rnd(gap_min_size, gap_max_size);
            for(ull i=0; i<gap_size; ++i) {
                output_file << 'N';
            }
        } // leading gap

        ull num_contigs = rnd(min_num_contigs, max_num_contigs);
        for(ull c=0; c<num_contigs; ++c) {
            ull contig_size = rnd(contig_min_size, contig_max_size);
            for(ull i=0; i<contig_size; ++i) {
                output_file << rndACGT();
            }

            if(c == num_contigs-1 && std::rand()%2 == 1) { // trailing gap
                ull gap_size = rnd(gap_min_size, gap_max_size);
                for(ull i=0; i<gap_size; ++i) {
                    output_file << 'N';
                }
            }
        }

        output_file << std::endl;
    }

    return EXIT_SUCCESS;
}