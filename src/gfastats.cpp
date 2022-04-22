//
//fastats.cpp
//
//Created by Giulio Formenti on 12/17/21.
//

#include <gfastats.h>

std::string version = "1.2.1";

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    unsigned long long int gSize = 0; // expected genome size, with 0 NG/LG* statistics are not computed
    int splitLength = 0; // line length for fasta output
    
    std::string iSeqFileArg; // input file to evaluate
    std::string iSakFileArg; // input of instructions for the swiss army knife
    std::string iAgpFileArg; // input agp
    std::string iBedIncludeFileArg; // input bed file of coordinates to include
    std::string iBedExcludeFileArg; // input bed file of coordinates to exclude
    
    std::string sortType = "none"; // type of sorting (default: none)
    
    std::string outSeq = "fasta"; // default output type
    
    char sizeOutType = 's'; // default output type with size flag (scaffold)
    char bedOutType = 'a'; // default output type with bed flag (agp)
    
    BedCoordinates bedInclude; // if include coordinates are provided as positional argument
    
    bool isPipe = false; // to check if input is from pipe
    char pipeType = 'n'; // default pipe type null
    
    if (argc == 1) { // gfastats with no arguments
            
        printf("gfastats input.[fasta|fastq|gfa][.gz] [expected genome size] [header[:start-end]]\n-h for additional help.\n");
        exit(0);
        
    }
    
    static struct option long_options[] = { // struct mapping long options
        {"input-sequence", required_argument, 0, 'f'},
        
        {"agp-to-path", required_argument, 0, 'a'}, // agp to path conversion
        {"swiss-army-knife", required_argument, 0, 'k'}, // the swiss army knife
        {"remove-terminal-gaps", no_argument, &rmGaps_flag, 1}, // this remove all gap edges at the end of sequences
        {"homopolymer-compress", required_argument, &hc_flag, 1},
        {"discover-paths", no_argument, &discoverPaths_flag, 1},
        {"sort", required_argument, 0, 0},
        
        {"include-bed", required_argument, 0, 'i'},
        {"exclude-bed", required_argument, 0, 'e'},
 
        {"out-format", required_argument, 0, 'o'},
        {"line-length", required_argument, 0, 0},
        {"out-sequence", no_argument, &outSequence_flag, 1},
        {"out-size", required_argument, 0, 's'},
        {"out-coord", required_argument, 0, 'b'},
        {"out-bubbles", no_argument, &outBubbles_flag, 1},
        
        {"stats", no_argument, &stats_flag, 1},
        {"seq-report", no_argument, &seqReport_flag, 1},
        {"nstar-report", no_argument, &nstarReport_flag, 1},
        {"tabular", no_argument, 0, 't'},
        {"locale", required_argument, 0, 0},
        
        {"verbose", no_argument, &verbose_flag, 1},
        {"cmd", no_argument, &cmd_flag, 1},
        {"version", no_argument, 0, 'v'},
        {"help", no_argument, 0, 'h'},
        
        {0, 0, 0, 0}
    };
    
    while (1) { // loop through argv
        
        int option_index = 0;
        
        c = getopt_long(argc, argv, "-:a:b:e:f:i:k:o:s:tvh",
                        long_options, &option_index);

        if (optind < argc && !isPipe) { // if pipe wasn't assigned already
            
            isPipe = isDash(argv[optind]) ? true : false; // check if the argument to the option is a '-' and set it as pipe input
            
        }
        
        if (optarg != nullptr && !isPipe) { // case where pipe input is given as positional argument (input sequence file)
        
            isPipe = isDash(optarg) ? true : false;
            
        }

        if (c == -1) { // exit the loop if run out of options
            break;
            
        }
        
        switch (c) {
            case ':': // handle options without arguments
                switch (optopt) { // the command line option last matched
                    case 'b':
                        bedOutType = 'a'; // default bed output is agp is -b option is given without argument
                        outCoord_flag = 1;
                        break;
                        
                    case 's':
                        sizeOutType = 's'; // default size output is scaffold is -s option is given without argument
                        outSize_flag = 1;
                        break;
                        
                    case 'o':
                        outSeq = "fasta"; // default output is fasta is -o option is given without argument
                        outFile_flag = 1;
                        break;
                        
                    default:
                        fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                        return EXIT_FAILURE;
                }
                break;
            default: // handle positional arguments
                if (pos_op == 1) { // only one positional argument given
                    
                    if (isPipe && pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                    
                        pipeType = 's'; // pipe input is a sequence
                    
                    }else{ // input is a regular file
                        
                        ifFileExists(optarg);
                        iSeqFileArg = optarg;
                        
                    }
                    
                    pos_op++;
                    
                }else if (pos_op == 2 || pos_op == 3) { // if >1 positional argument, check what additional positional arguments are present
                    
                    if (isInt(optarg)) { // if the positional argument is a number, it is likely the expected genome size
                        
                        gSize = atoll(optarg); pos_op++;
                        
                    }else{ // else it is an include argument
                        
                        std::string header = optarg, cBegin, cEnd; // the header for coordinates provided as positional argument
                        
                        reverse(header.begin(), header.end()); // we work our way from the end
                        
                        cBegin = header.substr(header.find('-') + 1, header.find(':') - header.find('-') - 1);
                        cEnd = header.substr(0, header.find('-'));
                        
                        if(isNumber(cEnd) &&
                           isNumber(cBegin)) { // prevent headers with : - characters to break the extraction
                            
                            header = header.substr(header.find(':') + 1, header.size());
                            reverse(header.begin(), header.end());
                            
                            reverse(cBegin.begin(), cBegin.end());
                            reverse(cEnd.begin(), cEnd.end());
                            
                        }else{
                            
                            reverse(header.begin(), header.end());
                            cBegin = "0";
                            cEnd = "0";
                            
                        }
                        
                        bedInclude.pushCoordinates(header, stoi(cBegin), stoi(cEnd)); pos_op++;
                        
                    }
                    
                }
                else{printf("Error: too many positional arguments (%s).\n",optarg);exit(1);}
                
            case 0: // case for long options without short options
                
                if (strcmp(long_options[option_index].name,"line-length") == 0)
                  splitLength = atoi(optarg);
                
                if (strcmp(long_options[option_index].name,"sort") == 0) {

                    std::vector<std::string> options {"none", "ascending", "descending", "largest", "smallest"};
                            
                    if (std::find(options.begin(), options.end(), optarg) != options.end() || ifFileExists(optarg)){
                        
                        sortType = optarg;
                        
                    }else{printf("Error: unrecognized sorting option (%s).\n", optarg); exit(1);}
                    
                }

                if(strcmp(long_options[option_index].name,"homopolymer-compress") == 0) {
                    hc_cutoff = atoi(optarg);
                    stats_flag = true;
                }
                
                if (strcmp(long_options[option_index].name,"locale") == 0) {
                    
                    setlocale(LC_ALL, optarg);
                    std::cout.imbue(std::locale(optarg));
                    std::locale::global(std::locale(optarg));
                    stats_flag = true;
                    
                }
                
                break;
                
            case 'a': // agp to paths
                
                if (isPipe && pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    pipeType = 'a'; // pipe input is agp
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    iAgpFileArg = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
            
            case 'b': // output bed type (agp, contig, gaps)
                bedOutType = *optarg;
                outCoord_flag = 1;
                break;
                
            case 'e': // bed exclude
                
                if (isPipe && pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    pipeType = 'e'; // pipe input is an exclude bed
                
                }else{ // input is a regular file
                
                    ifFileExists(optarg);
                    iBedExcludeFileArg = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
                
            case 'f': // input sequence
                
                if (isPipe && pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    pipeType = 's'; // pipe input is a sequence
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    iSeqFileArg = optarg;
                    stats_flag = true;
                    
                }
                    
                break;
                
            case 'i': // bed include
                
                if (isPipe && pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    pipeType = 'i'; // pipe input is an include bed
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    iBedIncludeFileArg = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
                
            case 'k': // the swiss army knife
                
                if (isPipe && pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                
                    pipeType = 'k'; // pipe input is a set of instructions for the swiss army knife
                
                }else{ // input is a regular file
                    
                    ifFileExists(optarg);
                    iSakFileArg = optarg;
                    
                }
                    
                stats_flag = 1;
                break;
                
            case 'o': // handle output (file or stdout)
                outSeq = optarg;
                outFile_flag = 1;
                break;
                
            case 's': // output size of features
                sizeOutType = *optarg;
                outSize_flag = 1;
                break;
                
            case 't': // tabular output
                tabular_flag = 1;
                break;
                
            case 'v': // software version
                printf("gfastats v%s.\n", version.c_str());
                printf("Giulio Formenti giulio.formenti@gmail.com.\n");
                exit(0);
                
            case 'h': // help
                printf("gfastats input.[fasta|fastq|gfa][.gz] [expected genome size] [header[:start-end]]\n");
                printf("genome size: estimated genome size for NG* statistics (optional).\n");
                printf("header: target specific sequence by header, optionally with coordinates (optional).\n");
                printf("\nOptions:\n");
                printf("-a --agp-to-path <file> converts input agp to path and replaces existing paths.\n");
                printf("-b --out-coord a|s|c|g generates bed coordinates of given feature (agp|scaffolds|contigs|gaps default:agp).\n");
                printf("-e --exclude-bed <file> opposite of --include-bed. They can be combined (no coordinates).\n");
                printf("-f --fasta <file> input file (fasta, fastq, gfa [.gz]). Also as first positional argument.\n");
                printf("-h --help print help and exit.\n");
                printf("-i --include-bed <file> generates output on a subset list of headers or coordinates in 0-based bed format.\n");
                printf("-k --swiss-army-knife <file> set of instructions provided as an ordered list.\n");
                printf("-o --out-format fasta|fastq|gfa[.gz] outputs selected sequences. If more than the extension is provided the output is written to the specified file (e.g. out.fasta.gz).\n");
                printf("-s --out-size s|c|g  generates size list of given feature (scaffolds|contigs|gaps default:scaffolds).\n");
                printf("-t --tabular output in tabular format.\n");
                printf("-v --version software version.\n");
                printf("--cmd print $0 to stdout.\n");
                printf("--discover-paths prototype to induce paths from input.\n");
                printf("--homopolymer-compress <threshhold> compress all the homopolymers in the input above the given threshhold.\n");
                printf("--line-length <n> specifies line length in when output format is fasta. Default has no line breaks.\n");
                printf("--nstar-report generates full N* and L* statistics.\n");
                printf("--out-sequence reports also the actual sequence (in combination with --seq-report).\n");
                printf("--out-bubbles outputs a potential list of bubbles in the graph.\n");
                printf("--seq-report report statistics for each sequence.\n");
                printf("--sort ascending|descending|largest|smallest|file sort sequences according to input. Ascending/descending used the sequence/path header.\n");
                printf("--stats report summary statistics (default).\n");
                printf("--verbose verbose output.\n");
                printf("--locale set a different locale, for instance to use , for thousand separators use en_US.UTF-8.\n");
                printf("\nAll input files can be piped from stdin using '-'.\n");
                exit(0);
        }
        
        if    (argc == 2 || // handle various cases in which the output should include summary stats
              (argc == 3 && pos_op == 2) ||
              (argc == 4 && pos_op == 3) ||
               nstarReport_flag ||
               discoverPaths_flag) {
            
            stats_flag = 1; // default mode 'stats'
            
        }
        
    }
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
    InFile inFile; // initialize sequence input file object
    
    verbose("File object generated");
    
    InSequences inSequences; // initialize sequence collection object
    
    verbose("Sequence object generated");
    
    inSequences = inFile.readFiles(iSeqFileArg, iSakFileArg, iAgpFileArg, iBedIncludeFileArg, iBedExcludeFileArg, bedInclude, isPipe, pipeType, sortType); // read the sequence input file object into the sequence collection object
    
    verbose("Finished reading sequences from file to sequence object");
    
    InSegment inSegment; // initialize a single input sequence object for output purposes
    
    Report report;
    
    if (seqReport_flag || outSequence_flag) { // report results for each sequence
        
        stats_flag = false;
        
        report.seqReport(inSequences, inSegment, outSequence_flag);
        
    }
    
    if (outFile_flag) { // output sequences to file or stdout
        
        stats_flag = false;
        
        report.outFile(inSequences, inSegment, splitLength, outSeq);
        
    }
    
    if (outSize_flag) { // output sequence sizes
        
        stats_flag = false;
        
        report.outSize(inSequences, inSegment, sizeOutType);
        
    }
    
    if (outCoord_flag) { // output coordinates
        
        stats_flag = false;
        
        report.outCoord(inSequences, inSegment, bedOutType);
        
    }
    
    if (stats_flag) { // output summary statistics
        
        report.reportStats(inSequences, gSize, bedOutType);
        
    }
    
    if (nstarReport_flag) { // output full N/L* statistics
        
        report.nstarReport(inSequences, gSize);
        
    }
    
    verbose("Generated output");
    
    exit(EXIT_SUCCESS);
    
}
