#ifndef GFA_LINES_H
#define GFA_LINES_H

class InSegment { // DNA sequence with no gaps
private:
    std::string seqHeader;
    std::string seqComment;
    std::string* inSequence = NULL;
    std::string* inSequenceQuality = NULL;
    unsigned long long int A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
    unsigned int uId = 0, iId = 0, seqPos = 0;
    std::vector<Tag> tags;
    
    friend class SAK;
    friend class InSequences;
    friend class Report;
    
public:
    
    ~InSegment();
    
    void set(Log* threadLog, unsigned int uId, unsigned int iId, std::string seqHeader, std::string* seqComment, std::string* sequence, unsigned long long int* A, unsigned long long int* C, unsigned long long int* G, unsigned long long int* T, unsigned long long int* lowerCount, unsigned int seqPos, std::string* sequenceQuality = NULL, std::vector<Tag>* inSequenceTags = NULL);
    
    void setSeqHeader(std::string* h);
    
    void setSeqComment(std::string c);
    
    void setInSequence(std::string* s);
    
    void setInSequenceQuality(std::string* q);
    
    void setSeqTags(std::vector<Tag>* t);

    void setuId(unsigned int i); // absolute id
    
    void setiId(unsigned int i); // temporary id, internal to scaffold
    
    void setSeqPos(unsigned int i); // temporary id, internal to scaffold

    std::string getSeqHeader();
    
    std::string getSeqComment();
    
    std::vector<Tag> getTags();
    
    std::string getInSequence(unsigned int start = 0, unsigned int end = 0);
    
    std::string getInSequenceQuality(unsigned int start = 0, unsigned int end = 0);
    
    unsigned int getSeqPos();
    
    unsigned long long int getSegmentLen(unsigned long long int start = 0, unsigned long long int end = 0);
    
    unsigned int getuId(); // absolute id
    
    unsigned int getiId(); // temporary id, internal to scaffold
    
    void setACGT(unsigned long long int* a, unsigned long long int* c, unsigned long long int* g, unsigned long long int* t);
    
    void setLowerCount(unsigned long long int* C);
    
    unsigned long long int getA();
    
    unsigned long long int getC();
    
    unsigned long long int getG();
    
    unsigned long long int getT();
    
    unsigned int getLowerCount(unsigned long long int start = 0, unsigned long long int end = 0);
    
    double computeGCcontent();
    
    bool trimSegment(unsigned int start, unsigned int end);
    
    bool rvcpSegment();
    
    bool invertSegment();
    
};

class InGap {
private:
//    unsigned long long int lineN; // useful if we wish to sort as is the original input
    std::string gHeader;
    char sId1Or, sId2Or;
    unsigned int uId, iId, sId1, sId2, dist;
    std::vector<Tag> tags;
    
    friend class SAK;
    friend class InSequences;
    
public:
    void newGap(unsigned int uId, unsigned int sId1, unsigned int sId2, const char& sId1or, const char& sId2or, unsigned int& dist, std::string gHeader = "", std::vector<Tag> tags = {});

    void setuId(unsigned int i); // absolute id
    
    void setiId(unsigned int i); // temporary id, internal to scaffold

    void setsId1(unsigned int i);
    
    void setsId2(unsigned int i);
    
    void setDist(unsigned int i);
    
    std::string getgHeader();
    
    unsigned int getuId();
    
    unsigned int getsId1();
    
    char getsId1Or();
    
    unsigned int getsId2();
    
    char getsId2Or();
    
    unsigned int getDist(unsigned int start = 0, unsigned int end = 0);
    
    std::vector<Tag> getTags();
    
    
    
};
class InEdge {
    private:
//    unsigned long long int lineN; // useful if we wish to sort as is the original input
    std::string cigar, eHeader;
    char sId1Or, sId2Or;
    unsigned int eUId, eId, sId1, sId2;
    std::vector<Tag> tags;
    
    friend class SAK;
    friend class InSequences;
    
public:
    void newEdge(unsigned int eUId, unsigned int sId1, unsigned int sId2, const char& sId1Or, const char& sId2Or, std::string cigar = "", std::string eHeader = "", std::vector<Tag> tags = {});
    
    bool operator==(const InEdge& e) const;

    void seteUId(unsigned int i); // absolute id
    
    void seteId(unsigned int i); // temporary id, internal to scaffold

    void setsId1(unsigned int i);
    
    void setsId2(unsigned int i);
    
    void setSeqTags(std::vector<Tag>* t);
    
    std::string getCigar();
    
    unsigned int geteUId();

    unsigned int geteId();
    
    unsigned int getsId1();
    
    char getsId1Or();
    
    unsigned int getsId2();
    
    char getsId2Or();
    
    std::vector<Tag> getTags();
    
};

class InPath {
    
private:
//    unsigned long long int lineN; // useful if we wish to sort as is the original input
    std::string pHeader, pComment;
    std::vector<PathComponent> pathComponents;
    unsigned int pUId, contigN = 0, seqPos;
    
    unsigned long long int length = 0, segmentLength = 0, lowerCount = 0, A = 0, C = 0, G = 0, T = 0;
    
    friend class SAK;
    friend class InSequences;

public:
    
    void newPath(unsigned int pUid, std::string h, std::string c = "", unsigned int seqpos = 0);

    void setpUId(unsigned int pUid);
    
    void setHeader(std::string pheader);
    
    void setComment(std::string c);
    
    void add(PathType type, unsigned int UId, char sign = '+', unsigned long long int start = 0, unsigned long long int end = 0);
    
    void append(std::vector<PathComponent> components);
    
    void clearPath();
    
    void setComponents(std::vector<PathComponent> newComponents);
    
    std::vector<PathComponent> getComponents();
    
    std::vector<PathComponent>* getComponentsByRef();

    unsigned int getpUId();
    
    std::string getHeader();
    
    std::string getComment();
    
    unsigned int getSeqPos();
    
    unsigned int getContigN();
    
    unsigned long long int getLen();
    
    unsigned long long int getA();
    
    unsigned long long int getC();
    
    unsigned long long int getG();
    
    unsigned long long int getT();
    
    unsigned long long int getSegmentLen();
    
    unsigned long long int getLowerCount();
    
    void revCom();
    
    void increaseContigN();
    
    void increaseGapN();
    
    void increaseLen(unsigned long long int n);
    
    void increaseSegmentLen(unsigned long long int n);
    
    void increaseLowerCount(unsigned long long int n);
    
    void increaseA(unsigned long long int n);
    
    void increaseC(unsigned long long int n);
    
    void increaseG(unsigned long long int n);
    
    void increaseT(unsigned long long int n);
    
};

#endif /* GFA_LINES_H */
