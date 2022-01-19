//
//fastats.cpp
//
//Created by Giulio Formenti on 12/17/21.
//

//global
static int verbose_flag;
static int seqReport_flag;
static int outSequence_flag;
static int nstarReport_flag;

#include <gfastats.h>

int main(int argc, char **argv) {
    
    static int outSize_flag;
    static int outCoord_flag;
    static int outFile_flag;
    static int stats_flag;
    static int cmd_flag;
    
    short int c;
    short unsigned int arg_counter;
    short unsigned int pos_op = 1;
    unsigned long long int gSize = 0;
    int splitLength = 0;
    
    std::string iSeqFileArg;
    std::string iBedIncludeFileArg;
    std::string iBedExcludeFileArg;
    
    std::string outSeq = "fasta";
    
    char sizeOutType = 's';
    char bedOutType = 'a';
    
    BedCoordinates bedInclude;
    std::string header;
    unsigned int cBegin = 0, cEnd = 0;
    char* coord;
    
    bool isPipe = false;
    char pipeType = 'n';
    
    if (argc == 1) {
            
            printf("gfastats input.[fasta|fastq|gfa][.gz] [genome size] [header[:start-end]]\n-h for additional help.\n");
            exit(0);
        
    }
    
    static struct option long_options[] = {
        {"fasta", required_argument, 0, 'f'},
        {"include-bed", required_argument, 0, 'i'},
        {"exclude-bed", required_argument, 0, 'e'},
 
        {"out-format", required_argument, 0, 'o'},
        {"line-length", required_argument, 0, 0},
        {"out-sequence", no_argument, &outSequence_flag, 1},
        {"out-size", required_argument, 0, 's'},
        {"out-coord", required_argument, 0, 'b'},
        
        {"stats", no_argument, &stats_flag, 1},
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
        
        c = getopt_long(argc, argv, "-:b:e:f:i:o:s:tvh",
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
        
        switch (c) {
            case ':':
                switch (optopt) {
                    case 'b':
                        bedOutType = 'a';
                        outCoord_flag = 1;
                        break;
                        
                    case 's':
                        sizeOutType = 's';
                        outSize_flag = 1;
                        break;
                        
                    case 'o':
                        outSeq = "fasta";
                        outFile_flag = 1;
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
                        iSeqFileArg = optarg;
                        
                    }
                    
                    pos_op++;
                    
                }else if (pos_op == 2 || pos_op == 3) {
                    
                    if (isInt(optarg)) {
                        
                        gSize = atoi(optarg); pos_op++;
                        
                    }else{
                        
                        header = std::string(strtok(strdup(optarg),":"));
                        
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
                
                if (strcmp(long_options[option_index].name,"line-length") == 0)
                  splitLength = atoi(optarg);
                
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
                    iSeqFileArg = optarg;
                    stats_flag = true;
                    
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
                
            case 'o':
                outSeq = optarg;
                outFile_flag = 1;
                break;
                
            case 's':
                sizeOutType = *optarg;
                outSize_flag = 1;
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
                printf("-o --out-format fasta|fastq|gfa[.gz] outputs selected sequences. If more than the extension is provided the output is written to the specified file (e.g. out.fasta.gz).\n");
                printf("\t--line-length <n> specifies line length in when output format is fasta. Default has no line breaks.\n");
                
                printf("-s --out-size s|c|g  generates bed coordinates of given feature (scaffolds|contigs|gaps default:scaffolds).\n");
                printf("-b a|c|g --out-coord generates bed coordinates of given feature (agp|contigs|gaps default:agp).\n");
                printf("-i --include-bed <file> generates output on a subset list of headers or coordinates in 0-based bed format.\n");
                printf("-e --exclude-bed <file> opposite of --include-bed. They can be combined (no coordinates).\n");
                printf("-t --tabular output in tabular format.\n");
                printf("-v --verbose verbose output.\n");
                printf("-h --help print help and exit.\n");
                printf("--stats report summary statistics (default).\n");
                printf("--seq-report report statistics for each sequence.\n");
                printf("\t--out-sequence reports also the actual sequence (in combination with --seq-report).\n");
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
    
    inSequences = inFile.readFiles(iSeqFileArg, iBedIncludeFileArg, iBedExcludeFileArg, bedInclude, isPipe, pipeType);
    
    verbose(verbose_flag, "Finished reading sequences from file to sequence object");
    
    InSequence inSequence;
    
    Report report;
    
    if (seqReport_flag || outSequence_flag) {
        
        stats_flag = true;
        
        report.seqReport(inSequences, inSequence, outSequence_flag);
        
    }
    
    if (outFile_flag) {
        
        stats_flag = false;
        
        report.outFile(inSequences, inSequence, splitLength, outSeq);
        
    }
    
    if (outSize_flag) {
        
        stats_flag = false;
        
        report.outSize(inSequences, inSequence, sizeOutType);
        
    }
    
    if (outCoord_flag) {
        
        stats_flag = false;
        
        report.outCoord(inSequences, inSequence, bedOutType);
        
    }
    
    if (stats_flag) {
        
        report.reportStats(inSequences, gSize, bedOutType);
        
    }
    
    if (nstarReport_flag) {
        
        report.nstarReport(inSequences, gSize);
        
    }
    
    verbose(verbose_flag, "Generated output");
    
    exit(EXIT_SUCCESS);
    
}
