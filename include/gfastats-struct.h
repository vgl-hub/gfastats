#ifndef GFASTATS_STRUCT
#define GFASTATS_STRUCT

struct Sequence {
    
    std::string header, comment, sequence, sequenceQuality;
    unsigned int seqPos;
    
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
