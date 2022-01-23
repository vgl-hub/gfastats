//
//  gfastats-classes.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef gfastatscommons_h
#define gfastatscommons_h

//classes
class BedCoordinates { // generic representation of bed coordinates
private:
    std::vector<std::string> seqHeaders;
    std::vector<unsigned int> cBegin;
    std::vector<unsigned int> cEnd;
    
public:
    
    void pushCoordinates(std::string h, unsigned int b = 0, unsigned int e = 0) { // reading coordinates
        
        seqHeaders.push_back(h);
        cBegin.push_back(b);
        cEnd.push_back(e);
        
    }
    
    bool empty() {
        
        return (seqHeaders.size()==0) ? true : false; // check if no coordinates are present
        
    }
    
    std::vector<std::string> getSeqHeaders() { // get all the headers
        
        return seqHeaders;
        
    }
    
    std::string getSeqHeader(unsigned int pos) { // get a specific header
        
        return seqHeaders[pos];
        
    }
    
    unsigned int getcBegin(unsigned int pos) { // get a specific start coordinate
        
        return cBegin[pos];
        
    }
    
    unsigned int getcEnd(unsigned int pos) { // get a specific end coordinate
        
        return cEnd[pos];
        
    }
    
};

class InSequence { // generic representation of a DNA sequence
private:
    std::string seqHeader;
    std::string seqComment;
    std::string inSequence;
    std::string inSequenceQuality;
    std::vector<unsigned int> contigBoundaries; // we store the coordinates of contig boundaries for easy access
    std::vector<unsigned int> gapBoundaries; // we store the coordinates of gap boundaries for easy access
    unsigned int A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
    friend class SAK;
    
public:
    
    void traverseInSequence(std::string* s) { // traverse the sequence to measure sequence properties
        
        unsigned int pos = 0, A = 0, C = 0, G = 0, T = 0, lowerCount = 0;
        bool wasN = false, pushbackGap = false;
        std::vector<unsigned int> caseBoundaries;
        std::vector<unsigned int> gapBoundaries;
        gapBoundaries.reserve(200);
        
        for (char &base : *s) {
            
            if (islower(base)) {
                
                lowerCount++;
                
            }
            
            switch (base) {
                    
                case 'N':
                case 'n':
                case 'X':
                case 'x': {
                    
                    if (!wasN) { // gap start
                        
                        pushbackGap = true;
                        
                    }
                    
                    if(pos == (s->length()-1)) { // end of scaffold
                        
                        if (!wasN){
                            
                            gapBoundaries.push_back(pos);
                            
                        }
                        
                        pushbackGap = true;
                        pos++;
                        
                    }
                    
                    wasN = true;
                    
                    break;
                }
                default: {
                    
                    if (wasN) { // internal gap end
                        
                        pushbackGap = true;
                        
                    }
                    
                    switch (base) {
                        case 'A':
                        case 'a':{
                            
                            A++;
                            break;
                            
                        }
                        case 'C':
                        case 'c':{
                            
                            C++;
                            break;
                            
                        }
                        case 'G':
                        case 'g': {
                            
                            G++;
                            break;
                            
                        }
                        case 'T':
                        case 't': {
                            
                            T++;
                            break;
                            
                        }
                            
                    }
                    
                    wasN = false;
                    
                }
                    
            }
            
            if (pushbackGap) {
                
                gapBoundaries.push_back(pos);
                pushbackGap = false;
                
            }
            
            pos++;
            
        }
        
        setSeqGapBoundaries(gapBoundaries);
        setSeqContigBoundaries(gapBoundaries);
        setACGT(A, C, G, T);
        setLowerCount(lowerCount);
        
    }
    
    void setSeqHeader(std::string* h) {
        seqHeader = *h;
    }
    
    void setSeqComment(std::string c) {
        seqComment = c;
    }
    
    void setInSequence(std::string* s) {
        inSequence = *s;
    }
    
    void setInSequenceQuality(std::string* q) {
        inSequenceQuality = *q;
    }
    
    void setSeqContigBoundaries(std::vector<unsigned int> &gapBoundaries) { // set the contig boundaries
        
        std::vector<unsigned int> newSeqContigBoundaries;
        
        newSeqContigBoundaries.reserve(gapBoundaries.size() + 2);
        
        if (gapBoundaries.size() > 0) {
            
            newSeqContigBoundaries = gapBoundaries;
            
            if (gapBoundaries[0] != 0) {
                
                newSeqContigBoundaries.insert(newSeqContigBoundaries.begin(), 0);
                
            }else{
                
                newSeqContigBoundaries.erase(newSeqContigBoundaries.begin());
                
            }
            
            if (newSeqContigBoundaries[newSeqContigBoundaries.size()-1] != inSequence.size()) {
                
                newSeqContigBoundaries.insert(newSeqContigBoundaries.end(), inSequence.size());
                
            }else{
                
                newSeqContigBoundaries.pop_back();
                
            }
            
        }else{
            
            newSeqContigBoundaries = {0, (unsigned int) inSequence.size()};
            
        }
        
        contigBoundaries = newSeqContigBoundaries;
        
    }
    
    void setSeqGapBoundaries(std::vector<unsigned int> &g) { // set gap boundaries
        gapBoundaries = g;
    }
    
    std::string getSeqHeader() {
        return seqHeader;
    }
    
    std::string getSeqComment() {
        return seqComment;
    }
    
    std::string getInSequence() {
        return inSequence;
    }
    
    std::string getInSequenceQuality() {
        
        return inSequenceQuality;
    }
    
    unsigned int getSeqScaffLen() {
        return inSequence.size();
    }
    
    std::vector<unsigned int> getSeqContigBoundaries() {
        return contigBoundaries;
    }
    
    std::vector<unsigned int> getSeqContigLens() {
        return intervalSizes(contigBoundaries);
    }
    
    unsigned int getContigSum() {
        
        unsigned int contigSum = 0;
        
        for (unsigned int& g : intervalSizes(contigBoundaries))
            contigSum += g;
        
        return contigSum;
    }
    
    unsigned int getContigN() {
        
        return intervalSizes(contigBoundaries).size();
    }
    
    std::vector<unsigned int> getSeqGapBoundaries() {
        return gapBoundaries;
    }
    
    std::vector<unsigned int> getSeqGapLens() {
        
        return intervalSizes(gapBoundaries);
    }
    
    unsigned int getGapSum() {
        
        unsigned int gapSum = 0;
        
        for (auto& g : intervalSizes(gapBoundaries))
            gapSum += g;
        
        return gapSum;
    }
    
    unsigned int getGapN() {
        
        return intervalSizes(gapBoundaries).size();
    }
    
    void setACGT(unsigned int a, unsigned int c, unsigned int g, unsigned int t) {
        
        A = a;
        C = c;
        G = g;
        T = t;
        
    }
    
    void setLowerCount(unsigned int C) {
        
        lowerCount = C;
        
    }
    
    unsigned int getA() {
        
        return A;
    }
    
    unsigned int getC() {
        
        return C;
    }
    
    unsigned int getG() {
        
        return G;
    }
    
    unsigned int getT() {
        
        return T;
    }
    
    unsigned int getLowerCount() {
        
        return lowerCount;
    }
    
    double computeGCcontent() {
        
        double GCcontent = (double) (G + C) / (G + C + A + T) * 100;
        
        return GCcontent;
    }
    
};

class InGap {
private:
    unsigned long long int lineN;
    std::string gId, sId1, sId2;
    char sId1Or, sId2Or;
    unsigned int dist;
    friend class SAK;
    friend class InSequences;
    int counter = 0;
    
public:
    bool readLine(std::string* newLine, unsigned long long int* lN) {
    
        lineN = *lN;
        
        gId = strtok(strdup((*newLine).c_str()),"\t"); // strip G
        
        gId = strtok(NULL,"\t");
        
        sId1 = strtok(NULL,"\t");
        sId2 = strtok(NULL,"\t");
        
        sId1Or = sId1.back(); // get sequence orientation in the gap
        sId2Or = sId2.back();
        
        sId1.pop_back(); // remove sequence orientation in the gap
        sId2.pop_back();
        
        dist = std::stoi(strtok(NULL,"\t"));
        
        counter++;
        
        return true;
        
    }
    
    std::string getgIds() {
        
        return gId;
        
    }
    
    std::string getsIds1() {
        
        return sId1;
        
    }
    
    std::string getsIds2() {
        
        return sId2;
        
    }
    
    unsigned int getsDists() {
        
        return dist;
        
    }
    
};
class InEdges {};
class InOlines {};
class InUlines {};

class InSequences { //collection of InSequence objects and their summary statistics
    
private:
    std::vector<InSequence> newSeq = std::vector<InSequence>();
    
    std::vector<unsigned int> scaffLens;
    std::vector<unsigned int> contigLens;
    std::vector<unsigned int> gapLens;
    
    std::vector<unsigned int> scaffNstars   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffLstars   {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffNGstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> scaffLGstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<unsigned int> contigNstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigLstars  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigNGstars {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> contigLGstars {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    std::vector<unsigned int> gapNstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    std::vector<unsigned int> gapLstars     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    double scaffAuN = 0, scaffAuNG = 0, contigAuN = 0, contigAuNG = 0, gapAuN = 0;
    
    InSequence inSequence;
    
    unsigned long long int totScaffLen = 0;
    unsigned long long int totContigLen = 0;
    
    unsigned int
    totGapLen = 0,
    totGapN = 0;
    
    unsigned long int totA = 0;
    unsigned long int totC = 0;
    unsigned long int totG = 0;
    unsigned long int totT = 0;
    unsigned long int totLowerCount = 0;
    
    //gfa variables
    std::vector<InGap> inGaps;
    std::vector<std::vector<Tuple>> adjList;
    std::vector<std::vector<Tuple>> adjListRV;
    std::unordered_map<std::string, unsigned long long int> ump;
    std::unordered_map<int, bool> visited;
    bool backward = false;
    //InEdges inEdges;
    //InOlines inOlines;
    //InUlines inUlines;
    
    friend class SAK;
    
public:
    
    //gfa methods
    void insertHash(std::string seqHeader, unsigned long long int seqN) {
    
        ump.insert({seqHeader, seqN});
    
    }
        
    void buildGraph(std::vector<InGap> const& edges) // graph Constructor
    {
        
        adjList.resize(newSeq.size()); // resize the vertex vector to hold `n` elements
        adjListRV.resize(newSeq.size()); // resize the vertex vector to hold `n` elements
 
        
        for (auto &edge: edges) // add edges to the graph
        {
            
            adjList[ump.at(edge.sId1)].push_back(std::make_tuple(edge.sId1Or, ump.at(edge.sId2), edge.sId2Or, edge.dist)); // insert at gap start gap destination and weight (gap size)
            adjListRV[ump.at(edge.sId2)].push_back(std::make_tuple(edge.sId2Or, ump.at(edge.sId1), edge.sId1Or, edge.dist)); // undirected graph

        }
        
    }
    
    void DFS(int v, InSequence &inSequenceNew, std::string &inSequence)
    {
    
        visited[v] = true; // mark the current node as visited
        std::string inSequenceNext;
        
        if (adjList[v].size() > 0 && adjListRV[v].size() > 0 && !backward) { // if there is more than one connected component to the vertex
            
            inSequence = (std::get<0>(adjList.at(v).at(0)) == '+') ? inSequence : reverse(inSequence); // check if vertex should be in forward orientation, if not reverse-complement
            
            inSequenceNext = (std::get<0>(adjListRV.at(v).at(0)) == '+') ? newSeq[v].getInSequence() : reverse(newSeq[v].getInSequence()); // check if vertex should be in forward orientation, if not reverse-complement
            
            inSequence += inSequenceNext;
                
        }else if (adjListRV[v].size() > 0 && !(adjList[v].size() > 0)){ // this is the final vertex
            
            inSequenceNext = (std::get<0>(adjListRV.at(v).at(0)) == '+') ? newSeq[v].getInSequence() : reverse(newSeq[v].getInSequence());
            
            inSequence += inSequenceNext;
            
            backward = true; // reached the end
            
        }else if (adjListRV[v].size() > 0){ // this is an intermediate vertex, only walking back
            
            inSequenceNext = (std::get<0>(adjListRV.at(v).at(0)) == '+') ? newSeq[v].getInSequence() : reverse(newSeq[v].getInSequence());
            
            inSequence.insert(0, inSequenceNext);
        
        }else if(adjList[v].size() == 0 && adjListRV[v].size() == 0){ // disconnected component
            
            inSequence += newSeq[v].getInSequence();
            
        }else{ // this is the first vertex
                
            inSequenceNext = (std::get<0>(adjList.at(v).at(0)) == '+') ? newSeq[v].getInSequence() : reverse(newSeq[v].getInSequence());
            
            inSequence.insert(0, inSequenceNext);
            
        }
        
        for (Tuple i: adjList[v]) { // recur for all forward vertices adjacent to this vertex
            if (!visited[std::get<1>(i)]) {
                
                inSequence += std::string(std::get<3>(i), 'N'); // add gaps
                
                DFS(std::get<1>(i), inSequenceNew, inSequence); // recurse
                
            }
        }
        
        for (Tuple i: adjListRV[v]) { // recur for all backward vertices adjacent to this vertex
            if (!visited[std::get<1>(i)]) {

                inSequence.insert(0,std::string(std::get<3>(i), 'N')); // add gaps

                DFS(std::get<1>(i), inSequenceNew, inSequence); // recurse

            }
        }
        
        inSequenceNew.setInSequence(&inSequence);
        
    }
    
    std::vector<std::vector<Tuple>> getAdjList() {
    
        return adjList;
        
    }
    
    bool getVisited(unsigned long long int seqN) {
    
        return visited[seqN];
        
    }
        
    //end of gfa methods
    
    void appendSequence(std::string* seqHeader, std::string* seqComment, std::string* sequence, std::string* sequenceQuality = NULL) { // method to append a new sequence
        
        inSequence.setSeqHeader(seqHeader);
        
        if (seqComment != NULL) {
            
            inSequence.setSeqComment(*seqComment);
            
        }
        
        verbose(verbose_flag, "Header, comment, and sequence read");
        
        verbose(verbose_flag, "Processing scaffold: " + *seqHeader);
        
        inSequence.setInSequence(sequence);
        
        verbose(verbose_flag, "Sequence set");
        
        inSequence.traverseInSequence(sequence);
        
        verbose(verbose_flag, "Sequence traversed");
        
        if (sequenceQuality != NULL) {
            
            inSequence.setInSequenceQuality(sequenceQuality);
            
            verbose(verbose_flag, "Sequence quality set");
            
        }
        
        newSeq.push_back(inSequence);
        
        verbose(verbose_flag, "Sequence added to sequence vector");
        
        increaseTotScaffLen(inSequence.getSeqScaffLen());
        
        verbose(verbose_flag, "Increased total scaffold length");
        
        recordScaffLen(inSequence.getSeqScaffLen());
        
        verbose(verbose_flag, "Recorded length of sequence");
        
        increaseTotContigLen(inSequence.getContigSum());
        
        verbose(verbose_flag, "Increased total contig length");
        
        recordContigLens(inSequence.getSeqContigLens());
        
        verbose(verbose_flag, "Recorded length of contigs in sequence");
        
        recordGapLens(inSequence.getSeqGapLens());
        
        verbose(verbose_flag, "Recorded length of gaps in sequence");
        
        increaseTotGapLen(inSequence.getGapSum());
        
        verbose(verbose_flag, "Increased total gap length");
        
        increaseGapN(inSequence.getGapN());
        
        verbose(verbose_flag, "Increased total number of gaps");
        
        increaseTotACGT(inSequence.getA(), inSequence.getC(), inSequence.getG(), inSequence.getT());
        
        verbose(verbose_flag, "Increased ACGT counts");
        
        increaseTotLowerCount(inSequence.getLowerCount());
        
        verbose(verbose_flag, "Increased count of upper bases");
        
        if(verbose_flag) {std::cout<<"\n";};
        
    }
    
    InSequence getInSequence(unsigned int &idx) {
        
        InSequence inSequence = newSeq[idx];
        return inSequence;
        
    }
    
    std::vector<InSequence> getInSequences() {
        
        return newSeq;
        
    }
    
    void increaseTotScaffLen(unsigned int ScaffLen) {
        
        totScaffLen += ScaffLen;
        
    }
    
    unsigned long long int getTotScaffLen() {
        
        return totScaffLen;
        
    }
    
    void increaseTotContigLen(unsigned int contigLen) {
        
        totContigLen += contigLen;
        
    }
    
    unsigned long long int getTotContigLen() {
        
        return totContigLen;
        
    }
    
    void increaseTotGapLen(unsigned int gapLen) {
        
        totGapLen += gapLen;
        
    }
    
    unsigned int getTotGapLen() {
        
        return totGapLen;
        
    }
    
    void increaseGapN(unsigned int gapN) {
        
        totGapN += gapN;
        
    }
    
    unsigned int getTotGapN() {
        
        return totGapN;
        
    }
    
    void recordScaffLen(unsigned int seqLen) {
        
        scaffLens.push_back(seqLen);
        
    }
    
    void recordContigLens(std::vector <unsigned int> seqLens) {
        
        std::vector <unsigned int> newContigLens;
        
        newContigLens.reserve(contigLens.size() + seqLens.size());
        newContigLens.insert(newContigLens.end(), contigLens.begin(), contigLens.end());
        newContigLens.insert(newContigLens.end(), seqLens.begin(), seqLens.end());
        
        contigLens = newContigLens;
        
    }
    
    void recordGapLens(std::vector <unsigned int> seqLens) {
        
        std::vector <unsigned int> newGapLens;
        
        newGapLens.reserve(gapLens.size() + seqLens.size());
        newGapLens.insert(newGapLens.end(), gapLens.begin(), gapLens.end());
        newGapLens.insert(newGapLens.end(), seqLens.begin(), seqLens.end());
        
        gapLens = newGapLens;
        
    }
    
    void evalNstars(char type, unsigned long long int gSize = 0) { // switch between scaffold, contig, gap while computing N* statistics
        
        switch(type) {
                
            case 's': {
                
                computeNstars(scaffLens, scaffNstars, scaffLstars, &scaffNGstars, &scaffLGstars, gSize);
                break;
                
            }
                
            case 'c': {
                
                computeNstars(contigLens, contigNstars, contigLstars, &contigNGstars, &contigLGstars, gSize);
                break;
                
            }
                
            case 'g': {
                
                computeNstars(gapLens, gapNstars, gapLstars);
                break;
                
            }
                
        }
        
    }
    
    
    void computeNstars(std::vector<unsigned int>& lens, // compute N/L* statistics, vector of all lengths
                       std::vector<unsigned int>& Nstars,      std::vector<unsigned int>& Lstars, // required arguments are passed by reference
                       std::vector<unsigned int>* NGstars = 0, std::vector<unsigned int>* LGstars = 0, unsigned long long int gSize = 0) { // optional arguments are passed by pointer
        
        sort(lens.begin(), lens.end(), std::greater<unsigned int>()); // sort lengths Z-A
        
        unsigned long long int sum = 0, totLen = 0;
        
        for(std::vector<unsigned int>::iterator it = lens.begin(); it != lens.end(); ++it) // find total length
            totLen += *it;
        
        short int N = 1, NG = 1;
        
        for(unsigned int i = 0; i < lens.size(); i++) { // for each length
            
            sum += lens[i]; // increase sum
            
            while (sum >= ((double) totLen / 10 * N) && N<= 10) { // conditionally add length.at or pos to each N/L* bin
                
                Nstars[N-1] = lens[i];
                Lstars[N-1] = i + 1;
                
                N = N + 1;
                
            }
                
            while (gSize > 0 && (sum >= ((double) gSize / 10 * NG)) && NG<= 10) { // if not computing gap statistics repeat also for NG/LG* statistics
                
                (*NGstars)[NG-1] = lens[i];
                (*LGstars)[NG-1] = i + 1;
                
                NG = NG + 1;
                
            }
            
        }
        
    }
    
    void evalAuN(char type, unsigned long long int gSize = 0) { // switch between scaffold, contig, gap while computing N* statistics
        
        switch(type) {
                
            case 's': {
                
                computeAuN(scaffLens, scaffAuN, &scaffAuNG, gSize);
                break;
                
            }
                
            case 'c': {
                
                computeAuN(contigLens, contigAuN, &contigAuNG, gSize);
                break;
                
            }
                
            case 'g': {
                
                computeAuN(gapLens, gapAuN);
                break;
                
            }
                
        }
        
    }
    
    void computeAuN(std::vector<unsigned int>& lens, double& auN, double* auNG = 0, unsigned long long int gSize = 0) {// compute N* statistics
        
        unsigned long long int totLen = 0;
        
        for(std::vector<unsigned int>::iterator it = lens.begin(); it != lens.end(); ++it) // find total length
            totLen += *it;
        
        for(unsigned int i = 0; i < lens.size(); i++) { // for each length
            
            auN += (double) lens[i] * lens[i] / totLen; // the area under the curve is the length (height) * fraction of the total length
            
            if (gSize > 0) {
                
                *auNG += (double) lens[i] * lens[i] / gSize;
                
            }
            
        }
        
    }
    
    unsigned int getScaffN() {
        
        return newSeq.size();
        
    }
    
    std::vector <unsigned int> getScaffNstars() {
        
        return scaffNstars;
        
    }
    
    std::vector <unsigned int> getScaffNGstars() {
        
        return scaffNGstars;
        
    }
    
    std::vector <unsigned int> getScaffLstars() {
        
        return scaffLstars;
        
    }
    
    std::vector <unsigned int> getScaffLGstars() {
        
        return scaffLGstars;
        
    }
    
    std::vector <unsigned int> getContigNstars() {
        
        return contigNstars;
        
    }
    
    std::vector <unsigned int> getContigNGstars() {
        
        return contigNGstars;
        
    }
    
    std::vector <unsigned int> getContigLstars() {
        
        return contigLstars;
        
    }
    
    std::vector <unsigned int> getContigLGstars() {
        
        return contigLGstars;
        
    }
    
    std::vector <unsigned int> getGapNstars() {
        
        return gapNstars;
        
    }
    
    std::vector <unsigned int> getGapLstars() {
        
        return gapLstars;
        
    }
    
    unsigned int getScaffN50() {
        
        return scaffNstars[4];
        
    }
    
    unsigned int getScaffNG50() {
        
        return scaffNGstars[4];
        
    }
    
    unsigned int getScaffL50() {
        
        return scaffLstars[4];
        
    }
    
    unsigned int getScaffLG50() {
        
        return scaffLGstars[4];
        
    }
    
    double getScaffauN() {
        
        return scaffAuN;
        
    }
    
    double getScaffauNG() {
        
        return scaffAuNG;
        
    }
    
    double getContigauN() {
        
        return contigAuN;
        
    }
    
    double getContigauNG() {
        
        return contigAuNG;
        
    }
    
    double getGapauN() {
        
        return gapAuN;
        
    }
    
    unsigned int getContigN() {
        
        return contigLens.size();
        
    }
    
    unsigned int getContigN50() {
        
        return contigNstars[4]; // middle value
        
    }
    
    unsigned int getContigNG50() {
        
        return contigNGstars[4];
        
    }
    
    unsigned int getContigL50() {
        
        return contigLstars[4];
        
    }
    
    unsigned int getContigLG50() {
        
        return contigLGstars[4];
        
    }
    
    unsigned int getGapN50() {
        
        return gapNstars[4];
        
    }
    
    unsigned int getGapL50() {
        
        return gapLstars[4];
        
    }
    
    unsigned int getLargestScaffold() {
            
        return scaffLens[0]; // sorted during N/L* computation
        
    }
    
    unsigned int getLargestContig() {
            
        return contigLens[0]; // sorted during N/L* computation
        
    }
    
    unsigned int getLargestGap() {
            
        return gapLens[0]; // sorted during N/L* computation
        
    }
    
    double computeAverageScaffLen() {
        
        return (double) totScaffLen/scaffLens.size();
        
    }
    
    double computeAverageContigLen() {
        
        return (double) totContigLen/contigLens.size();
        
    }
    
    double computeAverageGapLen() {
        
        return (double) totGapLen/gapLens.size();
        
    }
    
    void increaseTotACGT(unsigned int A, unsigned int C, unsigned int G, unsigned int T) {
        
        totA += A;
        totC += C;
        totG += G;
        totT += T;
        
    }
    
    void increaseTotLowerCount(unsigned int C) {
        
        totLowerCount += C;
        
    }
    
    unsigned long int getTotA() {
        
        return totA;
    }
    
    unsigned long int getTotC() {
        
        return totC;
    }
    
    unsigned long int getTotG() {
        
        return totG;
    }
    
    unsigned long int getTotT() {
        
        return totT;
    }
    
    unsigned long int getTotLowerCount() {
        
        return totLowerCount;
    }
    
    double computeGCcontent() {
        
        double GCcontent;
        
        if (newSeq.size()>0) {
            
            GCcontent = (double) (totC + totG) / (totA + totC + totG + totT) * 100;
            
        }else{
            
            GCcontent = 0;
            
        }
        
        return GCcontent;
    }
    
    //gfa methods
    bool appendGFAGap(InGap inGap) {
        
        inGaps.push_back(inGap);
        
        return true;
        
    }
    
    std::vector<InGap> getGFAGaps() { // return gfa gaps vector
        
        return inGaps;
        
    }
    
};

class SAK { // the swiss army knife
private:
    InSequences inSequences;
    InSequence inSequence1, inSequence2, inSequenceNew;
    std::string sId1Header, sId2Header;
    bool seq1 = false, seq2 = false;
    
public:
    
    bool joinByGap(InSequences &inSequences) { // joins two sequences via a gap based on instruction in gfa format
        
        std::vector<InGap> inGaps = inSequences.inGaps; // get the gaps vector
        
        for (std::vector<InGap>::const_iterator gapIt = inGaps.cbegin(); gapIt != inGaps.cend();) { // for each gap in the gap vector
            
            for (std::vector<InSequence>::const_iterator seqIt = inSequences.newSeq.cbegin(); seqIt != inSequences.newSeq.cend();) { //
                
                std::cout<<"sId1Header: "<<sId1Header<<std::endl;
                std::cout<<"sId2Header: "<<sId2Header<<std::endl;
                std::cout<<"(*seqIt).seqHeader: "<<(*seqIt).seqHeader<<std::endl;
                
                if ((*seqIt).seqHeader == sId1Header) {
                    
                    inSequence1 = (*seqIt);
                    seq1 = true;
                    inSequences.newSeq.erase(seqIt);
                    seqIt--;
                    
                }
                
                if ((*seqIt).seqHeader == sId2Header) {

                    inSequence2 = (*seqIt);
                    seq2 = true;
                    inSequences.newSeq.erase(seqIt);
                    seqIt--;
                    
                }
                
                if (seq1 && seq2) {break;}
                
                seqIt++;
                
            }
            
            seq1 = false, seq2 = false;
            
            inSequenceNew.seqHeader = (*gapIt).gId;
            inSequenceNew.seqComment = "JOIN " + inSequence1.seqHeader + " & " + inSequence2.seqHeader;
            inSequenceNew.inSequence = inSequence1.inSequence + std::string((*gapIt).dist, 'N') + inSequence2.inSequence;
            //inSequenceNew.inSequenceQuality = inSequence1.inSequenceQuality + std::string((*gapIt).dist, '!') + inSequence1.inSequenceQuality;

            inSequences.newSeq.push_back(inSequenceNew);
            
            gapIt++;
                
        }
        
        return true;
        
    }
    
};

#endif /* gfastats-commons_h */
