#ifndef GFASTATS_STRUCT
#define GFASTATS_STRUCT

struct UserInput { // a container for user input
  
    //files
    std::string iSeqFileArg; // input file to evaluate
    std::string iSakFileArg; // input of instructions for the swiss army knife
    std::string iAgpFileArg; // input agp
    std::string iBedIncludeFileArg; // input bed file of coordinates to include
    std::string iBedExcludeFileArg; // input bed file of coordinates to exclude
    std::string iReadFileArg; // input reads to evaluate
    
    //coordinates
    BedCoordinates bedIncludeList;
    
    //options
    char pipeType = 'n'; // default pipe type null
    std::string sortType = "none"; // type of sorting (default: none)
    
    std::string file(char type);
    
};

struct Sequence { // a generic sequence container
    
    std::string header, comment;
    std::string* sequence = NULL, *sequenceQuality = NULL;
    unsigned int seqPos = 0;
    
    ~Sequence();
    
};

struct Sequences { // a collection of sequences
    
    std::vector<Sequence*> sequences;
    unsigned int batchN;
    
    ~Sequences();
    
};

struct Tag {
    
    char type, label[3] = "";
    std::string content;
    
};

struct Gap {
    
    char orientation0;
    unsigned int segmentId;
    char orientation1;
    unsigned int dist;
    unsigned int edgeId;
    
};

struct Edge {
    
    char orientation0;
    unsigned int id;
    char orientation1;
    
    bool operator==(const Edge& e) const;
    
};

enum PathType { SEGMENT, GAP };
struct PathComponent {
    
    PathType type;
    unsigned int id;
    char orientation;
    unsigned long long int start;
    unsigned long long int end;
    
};

struct Bubble {
    unsigned int id0, id1, id2, id3;
};

#endif //GFASTATS_STRUCT
