//
//fastats.cpp
//
//Created by Giulio Formenti on 12/17/21.
//

//global
static int verbose_flag;
static int seqReport_flag;
static int nstarReport_flag;

#include <gfastats.h>

int main(int argc, char **argv) {
    
    static int outputSequence_flag;
    static int outputFasta_flag;
    static int stats_flag;
    static int cmd_flag;
    
    short int c;
    short unsigned int arg_counter;
    short unsigned int pos_op = 1;
    unsigned long long int gSize = 0;
    short int splitLength = 0;
    
    std::string iFastaFileArg;
    std::string iHeaderListFileArg;
    std::string iHeaderExcludeListFileArg;
    
    if (argc == 1) {
        printf("in.fasta [genome size] [target header]\n-h for additional help.\n");
        exit(0);
        
    }
    
    static struct option long_options[] = {
        {"fasta", required_argument, 0, 'f'},
        {"include-list", required_argument, 0, 'i'},
        {"exclude-list", required_argument, 0, 'e'},
        
        {"output-sequence", no_argument, &outputSequence_flag, 1},
        {"output-fasta", optional_argument, &outputFasta_flag, 1},
        
        {"stats", no_argument, 0, 's'},
        {"seq-report", no_argument, &seqReport_flag, 1},
        {"nstar-report", no_argument, &nstarReport_flag, 1},
        {"tabular", no_argument, 0, 't'},
        
        {"verbose", no_argument, 0, 'v'},
        {"cmd", no_argument, &cmd_flag, 1},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    
    while (1) {
        
        int option_index = 0;
        

        
        c = getopt_long(argc, argv, "-f:si:e:tvh",
                        long_options, &option_index);
            
        if (c == -1) {
            break;
            
        }
        
        if (outputFasta_flag && optarg != nullptr) {
            
            splitLength = atoi(optarg);
            
        }
        
        switch (c) {
            default:
                    if (pos_op == 1) {iFastaFileArg = optarg; pos_op++;}
                    else if (pos_op == 2 || pos_op == 3) {
                        
                        if (isInt(optarg)) {
                        
                            gSize = atoi(optarg); pos_op++;
                        
                        }else{
                            
                            headerList.push_back(optarg); pos_op++;
                            
                        }
                        
                    }
                    else{printf("Error: too many positional arguments (%s).\n",optarg);exit(1);}
            case 0:
                break;
                
            case 'f':
                iFastaFileArg = optarg;
                break;
                
            case 's':
                stats_flag = 1;
                break;
                
            case 'i':
                iHeaderListFileArg = optarg;
                stats_flag = 1;
                break;
                
            case 'e':
                iHeaderExcludeListFileArg = optarg;
                stats_flag = 1;
                break;

            case 't':
                tabular_flag = 1;
                break;
                
            case 'v':
                verbose_flag = 1;
                break;
                
            case 'h':
                printf("fastats in.fasta [genome size] [target header]\n");
                printf("genome size: estimated genome size for NG* statistics (optional).\n");
                printf("target header: compute statistics on a single header in the input (optional).\n\n");
                printf("Options:\n");
                printf("-f --fasta <file> fasta input. Also as first positional argument.\n");
                printf("-s --stats report summary statistics (default).\n");
                printf("-i --include-list <file> generates output on a subset list of headers.\n");
                printf("-e --exclude-list <file> opposite of --include-list. They can be combined.\n");
                printf("-t --tabular output in tabular format.\n");
                printf("-v --verbose verbose output.\n");
                printf("-h --help print help and exit.\n");
                printf("--seq-report report statistics for each sequence.\n");
                printf("--output-sequence reports also the actual sequence (in combination with --seq-report).\n");
                printf("--output-fasta [line length] generates a fasta output of the selected sequences. Default has no line breaks.\n");
                printf("--nstar-report generates full N* and L* statistics.\n");
                printf("--cmd print $0 to stdout.\n");
                exit(0);
        }
        
        if ( argc == 2 ||
             (argc == 3 && pos_op == 2) ||
             (argc == 4 && pos_op == 3) ||
             nstarReport_flag) {
            
            stats_flag = 1; // default mode 'stats'
            
        }
        
    }
    
    if (cmd_flag) {
        
        arg_counter = -1;
        while (arg_counter++ < argc-1) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
    FastaFile iFile;
    
    verbose(verbose_flag, "File object generated");
    
    FastaSequences fastaSequences;
    
    verbose(verbose_flag, "Fasta sequence object generated");
    
    fastaSequences = iFile.Read(iFastaFileArg, iHeaderListFileArg, iHeaderExcludeListFileArg);
    
    verbose(verbose_flag, "Finished reading sequences from file to fasta sequence object");
    
    unsigned int counter = 0;
    FastaSequence fastaSequence;
    
    if (seqReport_flag) {
        
        while (counter < fastaSequences.getScaffN()) {
            
            fastaSequence = fastaSequences.getFastaSequences(counter);
            
            std::cout<<output("Seq:")<<counter+1<<std::endl;
            std::cout<<output("Header:")<<fastaSequence.getFastaHeader()<<std::endl;
            std::cout<<output("Comment:")<<fastaSequence.getFastaComment()<<std::endl;
            std::cout<<output("Total sequence length:")<<fastaSequence.getFastaScaffLens()<<std::endl;
            std::cout<<output("Total contig length:")<<fastaSequence.getContigSum()<<std::endl;
            std::cout<<output("Total gap length:")<<fastaSequence.getGapSum()<<std::endl;
            std::cout<<output("Number of Gaps:")<<fastaSequence.getGapN()<<std::endl;
            
            printf("%s%u, %u, %u, %u\n",output("Base composition (ACGT):").c_str(), fastaSequence.getA(),
                   fastaSequence.getC(),
                   fastaSequence.getG(),
                   fastaSequence.getT());
            printf("%s%.2f\n",output("GC content %:").c_str(), fastaSequence.computeGCcontent());
            
            
            if (outputSequence_flag) {
                
                std::cout<<output("Sequence:")<<fastaSequence.getFastaSequence()<<std::endl;
                
            }
            
            std::cout<<std::endl;
            counter++;
            
        }
        
        counter = 0;
        
        std::cout<<output("+++Summary+++")<<std::endl;
        
    }
    
    if (outputFasta_flag) {
        
        stats_flag = false;
        
        while (counter < fastaSequences.getScaffN()) {
            
            fastaSequence = fastaSequences.getFastaSequences(counter);
            
            std::cout<<">"<<fastaSequence.getFastaHeader()<<" "<<fastaSequence.getFastaComment()<<std::endl;
            
            if (splitLength != 0) {
                
                unsigned int pos = 0;
                std::string line;
                
                for (char& base : fastaSequence.getFastaSequence())
                {
                    
                    line += base;
                    
                    if (pos == splitLength) {
                        
                        std::cout<<line;
                        std::cout<<std::endl;
                        
                        line = "";
                        pos = 0;
                        
                    }
                    
                    pos++;
                    
                }

                if (fastaSequence.getFastaSequence().length() % splitLength != 0) {
                    
                    std::cout<<std::endl;
                    
                }
                
            }else{
                
                std::cout<<fastaSequence.getFastaSequence()<<std::endl;
                
            }
            
            counter++;
            
        }
        
        
    }
    
    if (stats_flag) {
        
        verbose(verbose_flag, "Computed scaffN50");
        
        std::cout<<output("N scaffolds:")<<fastaSequences.getScaffN()<<std::endl;
        std::cout<<output("Total length:")<<fastaSequences.getTotScaffLen()<<std::endl;
        printf("%s%.2f\n",output("Average scaffold length:").c_str(), fastaSequences.computeAverageScaffLen());
        fastaSequences.computeScaffNstars(gSize);
        std::cout<<output("Scaffold N50:")<<fastaSequences.getScaffN50()<<std::endl;
        std::cout<<output("Scaffold L50:")<<fastaSequences.getScaffL50()<<std::endl;
        
        if (gSize > 0) {
            
            std::cout<<output("Scaffold NG50:")<<fastaSequences.getScaffNG50()<<std::endl;
            std::cout<<output("Scaffold LG50:")<<fastaSequences.getScaffLG50()<<std::endl;
            
        }
        std::cout<<output("N contigs:")<<fastaSequences.getContigN()<<std::endl;
        fastaSequences.computeContigNstars(gSize);
        std::cout<<output("Contig N50:")<<fastaSequences.getContigN50()<<std::endl;
        std::cout<<output("Contig L50:")<<fastaSequences.getContigL50()<<std::endl;
        
        if (gSize > 0) {
            
            std::cout<<output("Contig NG50:")<<fastaSequences.getContigNG50()<<std::endl;
            std::cout<<output("Contig LG50:")<<fastaSequences.getContigLG50()<<std::endl;
            
        }
        std::cout<<output("Largest scaffold:")<<fastaSequences.getLargestScaffold()<<std::endl;
        
        fastaSequences.computeGapNstars(gSize);
        std::cout<<output("Total gap length:")<<fastaSequences.getTotGapLen()<<std::endl;
        std::cout<<output("Number of Gaps:")<<fastaSequences.getTotGapN()<<std::endl;
        
        printf("%s%lu, %lu, %lu, %lu\n",output("Base composition (ACGT):").c_str(), fastaSequences.getTotA(),
               fastaSequences.getTotC(),
               fastaSequences.getTotG(),
               fastaSequences.getTotT());
        printf("%s%.2f\n",output("GC content %:").c_str(), fastaSequences.computeGCcontent());
        
        counter = 0;
        
    }
    
    if (nstarReport_flag) {
        
        int pos = 1;
        std::vector <unsigned int> scaffNstars = fastaSequences.getScaffNstars();
        for (unsigned int val : scaffNstars) {
            std::cout<<"Scaffold N"<<pos*10<<": "<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> scaffLstars = fastaSequences.getScaffLstars();
        for (unsigned int val : scaffLstars) {
            std::cout<<"Scaffold L"<<pos*10<<": "<<val<<std::endl;
            pos++;
        }
        
        if (gSize > 0) {
            
            pos = 1;
            std::vector <unsigned int> scaffNGstars = fastaSequences.getScaffNGstars();
            for (unsigned int val : scaffNGstars) {
                std::cout<<"Scaffold NG"<<pos*10<<": "<<val<<std::endl;
                pos++;
            }
            
            pos = 1;
            std::vector <unsigned int> scaffLGstars = fastaSequences.getScaffLGstars();
            for (unsigned int val : scaffLGstars) {
                std::cout<<"Scaffold LG"<<pos*10<<": "<<val<<std::endl;
                pos++;
            }
            
        }
        
        pos = 1;
        std::vector <unsigned int> contigNstars = fastaSequences.getContigNstars();
        for (unsigned int val : contigNstars) {
            std::cout<<"Contig N"<<pos*10<<": "<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> contigLstars = fastaSequences.getContigLstars();
        for (unsigned int val : contigLstars) {
            std::cout<<"Contig L"<<pos*10<<": "<<val<<std::endl;
            pos++;
        }
        
        if (gSize > 0) {
            
            pos = 1;
            std::vector <unsigned int> contigNGstars = fastaSequences.getContigNGstars();
            for (unsigned int val : contigNGstars) {
                std::cout<<"Contig NG"<<pos*10<<": "<<val<<std::endl;
                pos++;
            }
            
            pos = 1;
            std::vector <unsigned int> contigLGstars = fastaSequences.getContigLGstars();
            for (unsigned int val : contigLGstars) {
                std::cout<<"Contig LG"<<pos*10<<": "<<val<<std::endl;
                pos++;
            }
            
        }
        
        pos = 1;
        std::vector <unsigned int> gapNstars = fastaSequences.getGapNstars();
        for (unsigned int val : gapNstars) {
            std::cout<<"Gap N"<<pos*10<<": "<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> gapLstars = fastaSequences.getGapLstars();
        for (unsigned int val : gapLstars) {
            std::cout<<"Gap L"<<pos*10<<": "<<val<<std::endl;
            pos++;
        }
        
    }
    
    verbose(verbose_flag, "Generated output");
    
    exit(EXIT_SUCCESS);
    
}
