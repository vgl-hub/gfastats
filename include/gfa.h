//
//  gfastats-gfa.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFA_H
#define GFA_H

class InSequences { //collection of InSegment and inGap objects and their summary statistics
    
private:
    
    ThreadPool<std::function<void()>> threadPool;
    std::vector<Log> logs;
    
    //gfa variables
    std::vector<InSegment*> inSegments;
    std::vector<InGap> inGaps;
    std::vector<InEdge> inEdges;
    std::vector<InPath> inPaths;
    std::vector<InSegment*> inReads;
    std::vector<std::vector<Gap>> adjListFW;
    std::vector<std::vector<Gap>> adjListBW;
    std::vector<std::vector<Edge>> adjEdgeList;
    phmap::flat_hash_map<std::string, unsigned int> headersToIds;
    phmap::flat_hash_map<unsigned int, std::string> idsToHeaders;
    phmap::flat_hash_map<int, bool> visited, deleted;
    bool backward = false;
    
    std::vector<unsigned long long int> scaffLens;
    std::vector<unsigned long long int> contigLens;
    std::vector<unsigned long long int> gapLens;
    
    std::vector<unsigned long long int> scaffNstars   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffLstars   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned long long int> scaffNGstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffLGstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<unsigned long long int> contigNstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigLstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned long long int> contigNGstars {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigLGstars {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<unsigned long long int> gapNstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> gapLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<unsigned long long int> readNstars    {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    double scaffAuN = 0, scaffAuNG = 0, contigAuN = 0, contigAuNG = 0, gapAuN = 0;
    
    InSegment inSegment;
    InGap gap;
    InEdge edge;
    InPath path;
    
    unsigned long long int
    totScaffLen = 0,
    totContigLen = 0,
    totSegmentLen = 0,
    totGapLen = 0,
    totA = 0,
    totC = 0,
    totG = 0,
    totT = 0,
    totLowerCount = 0;
    
    unsigned int
    scaffN = 0,
    contigN = 0,
    pathN = 0;

    //connectivity
    unsigned int deadEnds = 0;
    unsigned int disconnectedComponents = 0;
    unsigned long long int lengthDisconnectedComponents = 0;
    
    std::vector<Bubble> bubbles;
    
    friend class SAK;
    
public:
    
    ~InSequences();
    
    UIdGenerator uId; // unique numeric identifier for each feature

    void threadPoolInit(int threadN);
    
    void threadStart(std::function<void()> job);

    bool threadEmpty();
    
    unsigned int threadQueueSize();
    
    void threadsJoin();
    
    std::vector<Log> getLogs();

    InSegment* addSegment(Log* threadLog, unsigned int uId, unsigned int iId, std::string seqHeader, std::string* seqComment, std::string* sequence, unsigned long long int* A, unsigned long long int* C, unsigned long long int* G, unsigned long long int* T, unsigned long long int* lowerCount, unsigned int seqPos, std::string* sequenceQuality = NULL, std::vector<Tag>* inSequenceTags = NULL);
    
    InGap pushbackGap(Log* threadLog, InPath* path, std::string* seqHeader, unsigned int* iId, unsigned int* dist, char sign, unsigned int uId1, unsigned int uId2);
    
    InSegment* pushbackSegment(unsigned int currId, Log* threadLog, InPath* path, std::string* seqHeader, std::string* seqComment, std::string* sequence, unsigned int* iId, unsigned long long int* A, unsigned long long int* C, unsigned long long int* G, unsigned long long int* T, unsigned long long int* lowerCount, unsigned long long int sStart, unsigned long long int sEnd, std::string* sequenceQuality = NULL);
    
    void traverseInSequence(Sequence* sequence);
    
    void traverseInSegment(Sequence* sequence, std::vector<Tag> inSequenceTags);
    
    void traverseInRead(Sequence* sequence);
    
    void appendSequence(Sequence* sequence);
    
    void appendSegment(Sequence* sequence, std::vector<Tag> inSequenceTags);
    
    void appendRead(Sequence* sequence);
    
    InSegment *getInSegment(unsigned int sId);
    
    std::vector<InSegment*>* getInSegments();
    
    std::vector<InGap>* getInGaps();
    
    InPath getInPath(unsigned int pId);
    
    std::vector<InPath> getInPaths();
    
    unsigned long long int getTotScaffLen();
    
    unsigned long long int getTotSegmentLen();
    
    void changeTotGapLen(unsigned int gapLen);
    
    unsigned int getGapNScaffold();
    
    unsigned int getGapN();

    unsigned int getEdgeN();
    
    unsigned int getPathN();
    
    void recordScaffLen(unsigned long long int seqLen);
    
    void recordGapLen(unsigned int gapLen);
    
    void evalNstars(char type, unsigned long long int gSize = 0);
    
    void computeNstars(std::vector<unsigned long long int>& lens, // compute N/L* statistics, vector of all lengths
                       std::vector<unsigned long long int>& Nstars,      std::vector<unsigned int>& Lstars, // required arguments are passed by reference
                       std::vector<unsigned long long int>* NGstars = 0, std::vector<unsigned int>* LGstars = 0, unsigned long long int gSize = 0);
    
    void evalAuN(char type, unsigned long long int gSize = 0);
    
    void computeAuN(std::vector<unsigned long long int>& lens, double& auN, double* auNG = 0, unsigned long long int gSize = 0);
    
    unsigned int getScaffN();
    
    unsigned int getTotContigN();
    
    std::vector <unsigned long long int> getScaffNstars();
    
    std::vector <unsigned long long int> getScaffNGstars();
    
    std::vector <unsigned int> getScaffLstars();
    
    std::vector <unsigned int> getScaffLGstars();
    
    std::vector <unsigned long long int> getContigNstars();
    
    std::vector <unsigned long long int> getContigNGstars();
    
    std::vector <unsigned int> getContigLstars();
    
    std::vector <unsigned int> getContigLGstars();
    
    std::vector <unsigned long long int> getGapNstars();
    
    std::vector <unsigned int> getGapLstars();
    
    unsigned long long int getScaffN50();
    
    unsigned long long int getScaffNG50();
    
    unsigned int getScaffL50();
    
    unsigned int getScaffLG50();
    
    double getScaffauN();
    
    double getScaffauNG();
    
    double getContigauN();
    
    double getContigauNG();
    
    double getGapauN();
    
    unsigned int getSegmentN();
    
    unsigned int getContigN50();
    
    unsigned int getContigNG50();
    
    unsigned int getContigL50();
    
    unsigned int getContigLG50();
    
    unsigned int getGapN50();
    
    unsigned int getGapL50();
    
    unsigned long long int getLargestScaffold();

    unsigned long long int getSmallestScaffold();

    
    unsigned long long int getLargestContig();

    unsigned long long int getSmallestContig();

    
    unsigned int getLargestGap();

    unsigned int getSmallestGap();
    
    double computeAvgScaffLen();
    
    unsigned long long int getTotContigLen();

    double computeAvgContigLen();
    
    double computeAvgSegmentLen();
    
    unsigned long long int getTotGapLen();
    
    double computeAverageGapLen();
    
    unsigned long long int getTotA();
    
    unsigned long long int getTotC();
    
    unsigned long long int getTotG();
    
    unsigned long long int getTotT();
    
    unsigned long long int getTotLowerCount();
    
    double computeGCcontent();
    
    //gfa methods
    bool addGap(InGap inGap);
    
    bool addPath(InPath path);
    
    std::vector<InGap> getGaps();

    std::vector<InEdge> getEdges();
    
    bool appendEdge(InEdge edge);
    
    //sorting methods

    void sortSegmentsByOriginal();
    
    void sortPathsByOriginal();
    
    void sortPathsByNameAscending();
    
    void sortPathsByNameDescending();
    
    void sortPathsByList(std::vector<std::string> headerList);
    
    void sortPathsBySize(bool largest);
    
    //gfa methods
    void insertHash(const std::string &segHeader, unsigned int i);
    
    unsigned int getuId();
    
    phmap::flat_hash_map<std::string, unsigned int>* getHash1();

    phmap::flat_hash_map<unsigned int, std::string>* getHash2();
    
    void buildGraph(std::vector<InGap> const& gaps);

    void buildEdgeGraph(std::vector<InEdge> const& edges);

    void dfsEdges(unsigned int v, unsigned int* componentLength);
    
    void dfsScaffolds(unsigned int v, unsigned int* scaffSize, unsigned int* A, unsigned int* C, unsigned int* G, unsigned int* T, unsigned int* lowerCount); // Depth First Search to explore graph connectivity
    
    std::vector<std::vector<Gap>> getAdjListFW();
    
    std::vector<std::vector<Gap>> getAdjListBW();
    
    bool getVisited(unsigned int uId);
    
    bool getDeleted(unsigned int uId);
    
    bool updateStats();
    
    bool removeTerminalGaps();

    unsigned int getDeadEnds();

    unsigned int getDisconnectedComponents();
    
    unsigned int getLengthDisconnectedComponents();
    
    // instruction methods
    
    std::vector<InGap> getGap(std::string* contig1, std::string* contig2 = NULL);
    
    std::vector<unsigned int> removeGaps(std::string* contig1, std::string* contig2 = NULL);
    
    bool deleteSegment(std::string* contig1);
    
    void removePath(unsigned int pUId, bool all = false, bool silent = false);
    
    void removeGap(unsigned int gUId, bool silent = false);
    
    void removePathsFromSegment(unsigned int uId);
    
    void removePathComponents(unsigned int uId);
    
    void removeSegmentInPath(unsigned int suId, InGap gap);
    
    void joinPaths(std::string pHeader, unsigned int pUId1, unsigned int pUId2, std::string gHeader, unsigned int gUId, char pId1Or, char pId2Or, unsigned int dist, unsigned int start1, unsigned int end1, unsigned int start2, unsigned int end2);
    
    InPath joinPathsByComponent(std::string seqHeader, unsigned int uId1, unsigned int uId2, unsigned int uId3);
    
    void splitPath(unsigned int guId, std::string pHeader1, std::string pHeader2);
    
    void clearPaths();
    
    void clearGaps();
    
    void renamePath(unsigned int pUId, std::string pHeader, unsigned int* newpUId);
    
    void revComPath(unsigned int pUId);
    
    void trimPathByUId(unsigned int pUId, unsigned int start, unsigned int end);

    void trimPathByRef(std::vector<PathComponent>& pathComponents, unsigned int start, unsigned int end);

    void trimPath(std::vector<PathComponent>* pathComponents, unsigned int start, unsigned int end);
    
    void trimComponent(PathComponent& component, int start, int end);
    
    int getComponentSize(PathComponent& component, bool original);
    
    unsigned int pathLen(unsigned int pUId);
    
    void walkPath(InPath* path);
    
    void discoverPaths();
    
    void dfsPath(unsigned int v, InPath& newPath); // Depth First Search to build a new path given a vertex
    
    void findBubbles();
    
    std::vector<Bubble>* getBubbles();
    
    // end of gfa methods
    
    // read methods
    
    unsigned int getReadN();
    
    unsigned long long int getTotReadLen();
    
    double computeAvgReadLen();
    
    unsigned long long int getReadN50();
    
    // end of read methods
    
};

#endif /* GFA_H */
