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
    
    static int outCoord_flag;
    static int outSequence_flag;
    static int outFasta_flag;
    static int stats_flag;
    static int cmd_flag;
    
    short int c;
    short unsigned int arg_counter;
    short unsigned int pos_op = 1;
    unsigned long long int gSize = 0;
    short unsigned int splitLength = 0;
    
    std::string iFastaFileArg;
    std::string iBedIncludeFileArg;
    std::string iBedExcludeFileArg;
    
    char bedOutType = 'a';
    
    BedCoordinates bedInclude;
    std::string header;
    unsigned int cBegin = 0, cEnd = 0;
    char* coord;
    
    bool isPipe = false;
    char pipeType;
    pipeType = 'n';
    
    if (argc == 1) {
            
            printf("gfastats input.[fasta|fastq|gfa][.gz] [genome size] [header[:start-end]]\n-h for additional help.\n");
            exit(0);
        
    }
    
    static struct option long_options[] = {
        {"fasta", required_argument, 0, 'f'},
        {"include-bed", required_argument, 0, 'i'},
        {"exclude-bed", required_argument, 0, 'e'},
        
        {"out-sequence", no_argument, &outSequence_flag, 1},
        {"out-coord", optional_argument, 0, 'b'},
        {"out-fasta", optional_argument, &outFasta_flag, 1},
        
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
        
        c = getopt_long(argc, argv, "-:b:e:f:i:stvh",
                        long_options, &option_index);

        if (optind < argc && !isPipe) {
            
            isPipe = isDash(argv[optind]) ? true : false;
            
        }
        
        if (optarg != nullptr && !isPipe) {
        
            isPipe = isDash(optarg) ? true : false;
            
        }

        if (c == -1) {
            break;
            
        }
        
        if (outFasta_flag && optarg != nullptr) {
            
            splitLength = atoi(optarg);
            
        }
        
        switch (c) {
            case ':':
                switch (optopt) {
                    case 'b':
                        bedOutType = 'a';
                        outCoord_flag = 1;
                        break;
                        
                    default:
                        fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                        return EXIT_FAILURE;
                }
                break;
            default:
                if (pos_op == 1) {
                    
                    if (isPipe && pipeType == 'n') {
                    
                        pipeType = 'f';
                    
                    }else{
                        
                        ifFileExists(optarg);
                        iFastaFileArg = optarg;
                        
                    }
                    
                    pos_op++;
                    
                }else if (pos_op == 2 || pos_op == 3) {
                    
                    if (isInt(optarg)) {
                        
                        gSize = atoi(optarg); pos_op++;
                        
                    }else{
                        
                        header = std::string(strtok(strdup(optarg),""));
                        
                        coord = strtok(NULL,"-");
                        
                        if (coord != NULL) {
                            
                            cBegin = atoi(coord);
                            
                            coord = strtok(NULL,"-");
                            
                            if (coord != NULL) {
                                
                                cEnd = atoi(coord);
                                
                            }else{printf("Error: missing end coordinate (%s).\n", header.c_str()); exit(1);}
                            
                        }
                        
                        bedInclude.pushCoordinates(header, cBegin, cEnd); pos_op++;
                        
                    }
                    
                }
                else{printf("Error: too many positional arguments (%s).\n",optarg);exit(1);}
            case 0:
                break;
                
            case 'b':
                bedOutType = *optarg;
                outCoord_flag = 1;
                break;
                
            case 'e':
                
                if (isPipe && pipeType == 'n') {
                
                    pipeType = 'e';
                
                }else{
                
                    ifFileExists(optarg);
                    iBedExcludeFileArg = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
                
            case 'f':
                
                if (isPipe && pipeType == 'n') {
                
                    pipeType = 'f';
                
                }else{
                    
                    ifFileExists(optarg);
                    iFastaFileArg = optarg;
                    
                }
                    
                break;
                
            case 'i':
                
                if (isPipe && pipeType == 'n') {
                
                    pipeType = 'i';
                
                }else{
                    
                    ifFileExists(optarg);
                    iBedIncludeFileArg = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
                
            case 's':
                stats_flag = 1;
                break;
                
            case 't':
                tabular_flag = 1;
                break;
                
            case 'v':
                verbose_flag = 1;
                break;
                
            case 'h':
                printf("gfastats input.[fasta|fastq|gfa][.gz] [genome size] [header[:start-end]]\n");
                printf("genome size: estimated genome size for NG* statistics (optional).\n");
                printf("header: target specific sequence by header, optionally with coordinates (optional).\n");
                printf("\nOptions:\n");
                printf("-f --fasta <file> fasta input. Also as first positional argument.\n");
                printf("-s --stats report summary statistics (default).\n");
                printf("-b c|g|a --out-coord generates bed coordinates of given feature (contigs|gaps|agp default:agp).\n");
                printf("-i --include-bed <file> generates output on a subset list of headers or coordinates in 0-based bed format.\n");
                printf("-e --exclude-bed <file> opposite of --include-bed. They can be combined.\n");
                printf("-t --tabular output in tabular format.\n");
                printf("-v --verbose verbose output.\n");
                printf("-h --help print help and exit.\n");
                printf("--seq-report report statistics for each sequence.\n");
                printf("--out-sequence reports also the actual sequence (in combination with --seq-report).\n");
                printf("--out-fasta [line length] generates a fasta output of the selected sequences. Default has no line breaks.\n");
                printf("--nstar-report generates full N* and L* statistics.\n");
                printf("--cmd print $0 to stdout.\n");
                printf("\nAll input files can be piped from stdin using '-'.\n");
                exit(0);
        }
        
        if   (argc == 2 ||
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
    
    InFile inFile;
    
    verbose(verbose_flag, "File object generated");
    
    InSequences inSequences;
    
    verbose(verbose_flag, "Sequence object generated");
    
    inSequences = inFile.readFiles(iFastaFileArg, iBedIncludeFileArg, iBedExcludeFileArg, bedInclude, isPipe, pipeType);
    
    verbose(verbose_flag, "Finished reading sequences from file to fasta sequence object");
    
    unsigned int counter = 0;
    InSequence fastaSequence;
    
    if (seqReport_flag || outSequence_flag) {
        
        while (counter < inSequences.getScaffN()) {
            
            fastaSequence = inSequences.getInSequence(counter);
            
            std::cout<<output("Seq")<<counter+1<<std::endl;
            std::cout<<output("Header")<<fastaSequence.getFastaHeader()<<std::endl;
            std::cout<<output("Comment")<<fastaSequence.getFastaComment()<<std::endl;
            std::cout<<output("Total sequence length")<<fastaSequence.getFastaScaffLen()<<std::endl;
            std::cout<<output("Total contig length")<<fastaSequence.getContigSum()<<std::endl;
            std::cout<<output("Total gap length")<<fastaSequence.getGapSum()<<std::endl;
            std::cout<<output("Number of gaps")<<fastaSequence.getGapN()<<std::endl;
            
            printf("%s%u, %u, %u, %u\n",output("Base composition (ACGT)").c_str(), fastaSequence.getA(),
                   fastaSequence.getC(),
                   fastaSequence.getG(),
                   fastaSequence.getT());
            printf("%s%.2f\n",output("GC content %").c_str(), fastaSequence.computeGCcontent());
            
            
            if (outSequence_flag) {
                
                std::cout<<output("Sequence")<<fastaSequence.getInSequence()<<std::endl;
                
            }
            
            std::cout<<std::endl;
            counter++;
            
        }
        
        counter = 0;
        
        std::cout<<output("+++Summary+++")<<std::endl;
        
    }
    
    if (outFasta_flag) {
        
        stats_flag = false;
        
        while (counter < inSequences.getScaffN()) {
            
            fastaSequence = inSequences.getInSequence(counter);
            
            std::cout<<">"<<fastaSequence.getFastaHeader()<<" "<<fastaSequence.getFastaComment()<<std::endl;
            
            if (splitLength != 0) {
                
                unsigned int pos = 0;
                std::string line;
                
                for (char& base : fastaSequence.getInSequence())
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
                
                if (fastaSequence.getInSequence().length() % splitLength != 0) {
                    
                    std::cout<<std::endl;
                    
                }
                
            }else{
                
                std::cout<<fastaSequence.getInSequence()<<std::endl;
                
            }
            
            counter++;
            
        }
        
        
    }
    
    if (outCoord_flag) {
        
        stats_flag = false;
        counter = 0;
        
        std::string fastaHeader;
        std::vector<unsigned int> fastaBoundaries;
        
        switch (bedOutType) {
                
            case 'c': {
                
                while (counter < inSequences.getScaffN()) {
                    
                    fastaSequence = inSequences.getInSequence(counter);
                    
                    fastaHeader = fastaSequence.getFastaHeader();
                    
                    fastaBoundaries = fastaSequence.getFastaContigBoundaries();
                    
                    std::vector<unsigned int>::const_iterator end = fastaBoundaries.cend();
                    
                    for (std::vector<unsigned int>::const_iterator it = fastaBoundaries.cbegin(); it != end;) {
                        
                        std::cout<<fastaHeader<<"\t"<<*it<<"\t"<<*(it+1)<<std::endl;
                        
                        it = it + 2;
                        
                    }
                    
                    counter++;
                    
                }
                
                break;
                
            }
                
            case 'g': {
                
                while (counter < inSequences.getScaffN()) {
                    
                    fastaSequence = inSequences.getInSequence(counter);
                    
                    fastaHeader = fastaSequence.getFastaHeader();
                    
                    fastaBoundaries = fastaSequence.getFastaGapBoundaries();
                    
                    std::vector<unsigned int>::const_iterator end = fastaBoundaries.cend();
                    
                    for (std::vector<unsigned int>::const_iterator it = fastaBoundaries.cbegin(); it != end;) {
                        
                        std::cout<<fastaHeader<<"\t"<<*it<<"\t"<<*(it+1)<<std::endl;
                        
                        it = it + 2;
                        
                    }
                    
                    counter++;
                    
                }
                
                break;
                
            }
                
            default:
            case 'a': {
                
                unsigned int ctgN = 1, item = 1, len = 0;
                
                while (counter < inSequences.getScaffN()) {
                    
                    fastaSequence = inSequences.getInSequence(counter);
                    unsigned int fastaScaffLen = fastaSequence.getFastaScaffLen();
                    
                    fastaHeader = fastaSequence.getFastaHeader();
                    
                    fastaBoundaries = fastaSequence.getFastaContigBoundaries();
                    
                    std::vector<unsigned int>::const_iterator begin = fastaBoundaries.cbegin();
                    std::vector<unsigned int>::const_iterator end = fastaBoundaries.cend();
                    auto last = std::prev(end);
                    
                    if (*begin>0) {
                        
                        std::cout<<fastaHeader<<"\t"<<1<<"\t"<<*begin<<"\t"<<1<<"\t"<<"N"<<"\t"<<*begin<<"\tscaffold\tyes\t"<<std::endl;
                        
                        item++;
                        
                    }
                    
                    for (std::vector<unsigned int>::const_iterator it = fastaBoundaries.cbegin(); it != end;) {
                        
                        len = *(it+1) - *it;
                        
                        std::cout<<fastaHeader<<"\t"<<*it+1<<"\t"<<*(it+1)<<"\t"<<item<<"\t"<<"W"<<"\t"<<fastaHeader+"."<<ctgN<<"\t1\t"<<len<<"\t+"<<std::endl;
                        
                        item++;
                        
                        if (ctgN != fastaBoundaries.size()/2) {
                            
                            len = *(it+2) - *(it+1);
                            
                            std::cout<<fastaHeader<<"\t"<<*(it+1)+1<<"\t"<<*(it+2)<<"\t"<<item<<"\t"<<"N"<<"\t"<<len<<"\tscaffold\tyes\t"<<std::endl;
                            
                            item++;
                            
                        }
                        
                        if (ctgN == fastaBoundaries.size()/2 && fastaScaffLen > *last) {
                            
                            len = fastaScaffLen - *(it+1);
                            
                            std::cout<<fastaHeader<<"\t"<<*(it+1)+1<<"\t"<<fastaScaffLen<<"\t"<<item<<"\t"<<"N"<<"\t"<<len<<"\tscaffold\tyes\t"<<std::endl;
                            
                            item++;
                            
                        }
                        
                        ctgN++;
                        it = it + 2;
                        
                    }
                    
                    ctgN = 1;
                    item = 1;
                    counter++;
                    
                }
                
                break;
            }
                
        }
        
    }
    
    if (stats_flag) {
        
        verbose(verbose_flag, "Computed scaffN50");
        
        std::cout<<output("N scaffolds")<<inSequences.getScaffN()<<std::endl;
        std::cout<<output("Total scaffold length")<<inSequences.getTotScaffLen()<<std::endl;
        printf("%s%.2f\n",output("Average scaffold length").c_str(), inSequences.computeAverageScaffLen());
        inSequences.computeScaffNstars(gSize);
        std::cout<<output("Scaffold N50")<<inSequences.getScaffN50()<<std::endl;
        std::cout<<output("Scaffold L50")<<inSequences.getScaffL50()<<std::endl;
        
        if (gSize > 0) {
            
            std::cout<<output("Scaffold NG50")<<inSequences.getScaffNG50()<<std::endl;
            std::cout<<output("Scaffold LG50")<<inSequences.getScaffLG50()<<std::endl;
            
        }
        std::cout<<output("Largest scaffold")<<inSequences.getLargestScaffold()<<std::endl;
        
        std::cout<<output("N contigs")<<inSequences.getContigN()<<std::endl;
        std::cout<<output("Total contig length")<<inSequences.getTotContigLen()<<std::endl;
        inSequences.computeContigNstars(gSize);
        std::cout<<output("Contig N50")<<inSequences.getContigN50()<<std::endl;
        std::cout<<output("Contig L50")<<inSequences.getContigL50()<<std::endl;
        
        if (gSize > 0) {
            
            std::cout<<output("Contig NG50")<<inSequences.getContigNG50()<<std::endl;
            std::cout<<output("Contig LG50")<<inSequences.getContigLG50()<<std::endl;
            
        }
        
        inSequences.computeGapNstars(gSize);
        std::cout<<output("N of Gaps")<<inSequences.getTotGapN()<<std::endl;
        std::cout<<output("Total gap length")<<inSequences.getTotGapLen()<<std::endl;
        
        printf("%s%lu, %lu, %lu, %lu\n",output("Base composition (ACGT)").c_str(), inSequences.getTotA(),
               inSequences.getTotC(),
               inSequences.getTotG(),
               inSequences.getTotT());
        printf("%s%.2f\n",output("GC content %").c_str(), inSequences.computeGCcontent());
        
        counter = 0;
        
    }
    
    if (nstarReport_flag) {
        
        int pos = 1;
        std::vector <unsigned int> scaffNstars = inSequences.getScaffNstars();
        for (unsigned int val : scaffNstars) {
            std::cout<<output("Scaffold N"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> scaffLstars = inSequences.getScaffLstars();
        for (unsigned int val : scaffLstars) {
            std::cout<<output("Scaffold L"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        if (gSize > 0) {
            
            pos = 1;
            std::vector <unsigned int> scaffNGstars = inSequences.getScaffNGstars();
            for (unsigned int val : scaffNGstars) {
                std::cout<<output("Scaffold NG"+std::to_string(pos*10))<<val<<std::endl;
                pos++;
            }
            
            pos = 1;
            std::vector <unsigned int> scaffLGstars = inSequences.getScaffLGstars();
            for (unsigned int val : scaffLGstars) {
                std::cout<<output("Scaffold LG"+std::to_string(pos*10))<<val<<std::endl;
                pos++;
            }
            
        }
        
        pos = 1;
        std::vector <unsigned int> contigNstars = inSequences.getContigNstars();
        for (unsigned int val : contigNstars) {
            std::cout<<output("Contig N"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> contigLstars = inSequences.getContigLstars();
        for (unsigned int val : contigLstars) {
            std::cout<<output("Contig L"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        if (gSize > 0) {
            
            pos = 1;
            std::vector <unsigned int> contigNGstars = inSequences.getContigNGstars();
            for (unsigned int val : contigNGstars) {
                std::cout<<output("Contig NG"+std::to_string(pos*10))<<val<<std::endl;
                pos++;
            }
            
            pos = 1;
            std::vector <unsigned int> contigLGstars = inSequences.getContigLGstars();
            for (unsigned int val : contigLGstars) {
                std::cout<<output("Contig LG"+std::to_string(pos*10))<<val<<std::endl;
                pos++;
            }
            
        }
        
        pos = 1;
        std::vector <unsigned int> gapNstars = inSequences.getGapNstars();
        for (unsigned int val : gapNstars) {
            std::cout<<output("Gap N"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
        pos = 1;
        std::vector <unsigned int> gapLstars = inSequences.getGapLstars();
        for (unsigned int val : gapLstars) {
            std::cout<<output("Gap L"+std::to_string(pos*10))<<val<<std::endl;
            pos++;
        }
        
    }
    
    verbose(verbose_flag, "Generated output");
    
    exit(EXIT_SUCCESS);
    
}
