//
//  gfastats-sak.h
//
//
//  Created by Giulio Formenti on 1/30/22.
//

#ifndef GFASTATS_SAK_H
#define GFASTATS_SAK_H

struct Instruction {
    std::string action;
    std::string path1;
    std::string path2;
    std::string path3;
    std::string contig1;
    std::string contig2;
    std::string comment1;
    std::string comment2;
    std::string gHeader = "";
    
    int compressThreshhold;
    char pId1Or, pId2Or, sId1Or, sId2Or;
    unsigned int gUId = 0, dist = 0, start1 = 0, end1 = 0, start2 = 0, end2 = 0;
};

class SAK;

struct Funcs {
public:
    void (SAK::*read)(Instruction&, std::vector<std::string>&);
    bool (SAK::*exec)(InSequences&, Instruction&);
};

class SAK { // the swiss army knife
public:
    void read(Instruction &instruction, std::vector<std::string> &arguments) {
        return;
    }

    void readJoin(Instruction &instruction, std::vector<std::string> &arguments) {
        size_t pos = 0, pos2 = 0, pos3 = 0;
        
        instruction.pId1Or = arguments[1].back(); // get sequence orientation in the gap
        arguments[1].pop_back(); // remove sequence orientation in the gap
        
        if(arguments[1].back() == ')') {
            pos = arguments[1].find('(');
            pos2 = arguments[1].find(':');
            pos3 = arguments[1].find(')');
            
            instruction.path1 = arguments[1].substr(0, pos);
            instruction.start1 = stoi(arguments[1].substr(pos + 1, pos2 - pos - 1));
            instruction.end1 = stoi(arguments[1].substr(pos2 + 1, pos3 - pos2 - 1));
        }else{
            instruction.path1 = arguments[1];
        }
    
        instruction.pId2Or = arguments[2].back(); // get sequence orientation in the gap
        arguments[2].pop_back(); // remove sequence orientation in the gap
        
        if(arguments[2].back() == ')') {
            pos = arguments[2].find('(');
            pos2 = arguments[2].find(':');
            pos3 = arguments[2].find(')');
            
            instruction.path2 = arguments[2].substr(0, pos);
            instruction.start2 = stoi(arguments[2].substr(pos + 1, pos2 - pos - 1));
            instruction.end2 = stoi(arguments[2].substr(pos2 + 1, pos3 - pos2 - 1));
        }else{
            instruction.path2 = arguments[2];
        }
        
        instruction.dist = stoi(arguments[3]);
        instruction.gHeader = arguments[4];
        instruction.path3 = arguments[5];
        
        if (arguments[6] != "") {
            instruction.gUId = stoi(arguments[6]);
        }
    }

    bool join(InSequences& inSequences, Instruction &instruction) { // joins two sequences via a gap based on instruction
        inSequences.joinPaths(instruction.path3, inSequences.headersToIds[instruction.path1], inSequences.headersToIds[instruction.path2], instruction.gHeader, instruction.gUId, instruction.pId1Or, instruction.pId2Or, instruction.dist, instruction.start1, instruction.end1, instruction.start2, instruction.end2); // generate a new path by joining the paths that contain the two segments
        return true;
    }

    void readSplit(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
        instruction.contig2 = arguments[2];
        instruction.path1 = arguments[3];
        instruction.path2 = arguments[4];
        instruction.comment1 = arguments[5];
        instruction.comment2 = arguments[6];
    }

    bool split(InSequences& inSequences, Instruction &instruction) { // splits two sequences removing the gap in between based on instruction
        std::vector<unsigned int> guIds = inSequences.removeGaps(&instruction.contig1, &instruction.contig2); // remove the gap
        inSequences.splitPath(guIds[0], instruction.path1, instruction.path2); // generate two new the paths splitting the original path
        return true;
    }

    void readExcise(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];

        if (arguments[2] != "") {
            instruction.dist = stoi(arguments[2]);
        }else{
            instruction.dist = 0;
        }
        instruction.gHeader = arguments[3];
    }

    bool excise(InSequences& inSequences, Instruction &instruction) { // excises a sequence, removing also edges if present and optionally adding a gap
        std::vector<InGap> oldGaps = inSequences.getGap(&instruction.contig1); // get the gaps associated with contig1
        InGap gap;
        inSequences.uId++;

        gap.newGap(inSequences.uId, oldGaps[0].getsId1(), oldGaps[1].getsId2(), oldGaps[0].getsId1Or(), oldGaps[1].getsId2Or(), instruction.dist, instruction.gHeader); // define the new gap
        
        inSequences.insertHash1(instruction.gHeader, inSequences.uId); // header to hash table
        inSequences.insertHash2(inSequences.uId, instruction.gHeader); // header to hash table
        
        inSequences.addGap(gap); // introduce the new gap
        inSequences.removeGaps(&instruction.contig1);
        inSequences.removeSegmentInPath(inSequences.headersToIds[instruction.contig1], gap); // removes the segment from the path
        
        InPath path;
        path.setHeader(instruction.contig1);
        path.add('S', inSequences.headersToIds[instruction.contig1], '+');
        inSequences.addPath(path);
            
        return true;
    }

    void readRemove(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
    }

    bool remove(InSequences& inSequences, Instruction &instruction) { // removes a segment
        inSequences.removePathsFromSegment(inSequences.headersToIds[instruction.contig1]); // remove the paths involving contig1
        inSequences.removeSegment(&instruction.contig1); // remove the segment
        inSequences.updateGapLens();
        return true;
    }

    void readErase(Instruction &instruction, std::vector<std::string> &arguments) {  
        size_t pos1 = 0, pos2 = 0;
        pos1 = arguments[1].find(":");
        pos2 = arguments[1].find("-");
        instruction.contig1 = arguments[1].substr(0, pos1);
        instruction.start1 = stoi(arguments[1].substr(pos1+1, pos2));
        instruction.end1 = stoi(arguments[1].substr(pos2+1, arguments[1].size()+1));
    }

    bool erase(InSequences& inSequences, Instruction &instruction) { // erases a portion of sequence
        inSequences.inSegments[inSequences.headersToIds[instruction.contig1]].trimSegment(instruction.start1, instruction.end1); // trim segment
        inSequences.changeTotSegmentLen(instruction.start1-instruction.end1);
        return true;
    }

    void readRvcp(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
    }

    bool rvcp(InSequences& inSequences, Instruction &instruction) { // reverse complement sequence
        unsigned int uId = inSequences.headersToIds[instruction.contig1], sIdx = 0;
        auto sId = find_if(inSequences.inSegments.begin(), inSequences.inSegments.end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node uId, find it

        if (sId != inSequences.inSegments.end()) {sIdx = std::distance(inSequences.inSegments.begin(), sId);} // gives us the segment index
        
        inSequences.inSegments[sIdx].rvcpSegment(); // rvcp segment
        
        return true;
    }

    void readInvert(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.contig1 = arguments[1];
    }

    bool invert(InSequences& inSequences, Instruction &instruction) { // invert sequence
        unsigned int uId = inSequences.headersToIds[instruction.contig1], sIdx = 0;
        auto sId = find_if(inSequences.inSegments.begin(), inSequences.inSegments.end(), [uId](InSegment& obj) {return obj.getuId() == uId;}); // given a node uId, find it

        if (sId != inSequences.inSegments.end()) {sIdx = std::distance(inSequences.inSegments.begin(), sId);} // gives us the segment index
        inSequences.inSegments[sIdx].invertSegment(); // invert segment
        
        return true;
    }

    void readCompress(Instruction &instruction, std::vector<std::string> &arguments) {
        instruction.compressThreshhold = stoi(arguments[0]);
    }

    bool compress(InSequences &inSequences, Instruction &instruction) {
        // auto inSegments = inSequences.getInSegments();
        
        // for(int i=0; i<inSegments->size(); ++i) {
            
        // }

        return true;
    }

    void readDecompress(Instruction &instruction, std::vector<std::string> &arguments) {

    }

    bool decompress(InSequences &inSequences, Instruction &instruction) {
        
        return true;
    }

private:
    InSequences inSequences;
    InSegment inSegment1, inSegment2, inSegmentNew;
    std::string sId1Header, sId2Header;

    // unordered map to handle out correspondence in following switch statement
    const phmap::flat_hash_map<std::string, Funcs> string_to_funcs { // std::pair<void (SAK::*)(Instruction&, std::vector<std::string>&), bool (SAK::*)(InSequences&, Instruction&)>
        {"JOIN", {&readJoin, &join} },
        {"SPLIT", {&readSplit, &split} },
        {"EXCISE", {&readExcise, &excise} },
        {"REMOVE", {&readRemove, &remove} },
        {"ERASE", {&readErase, &erase} },
        {"RVCP", {&readRvcp, &rvcp} },
        {"INVERT", {&readInvert, &invert} },
        {"COMPRESS", {&readCompress, &compress} },
        {"DECOMPRESS", {&readDecompress, &decompress} }
    };

public:
    
    Instruction readInstruction(std::string line) {
        
        std::string delimiter = "\t";
        std::vector<std::string> arguments;
        
        size_t pos = 0;
        
        while ((pos = line.find(delimiter)) != std::string::npos) {
            
            arguments.push_back(line.substr(0, pos));
            
            line.erase(0, pos + delimiter.length());
        
        }
        
        arguments.push_back(line); // last column
        
        Instruction instruction;
        
        for (auto & c: arguments[0]) instruction.action += (char) toupper(c);
        
        if(!string_to_funcs.count(instruction.action)) {
            fprintf(stderr, "unrecognized action %s\n", instruction.action.c_str());
            exit(EXIT_FAILURE);
        }

        (this->*(string_to_funcs.at(instruction.action).read))(instruction, arguments);
        
        verbose("Instruction read");
        
        return instruction;
        
    }
    
    bool executeInstruction(InSequences& inSequences, Instruction instruction) {
        if(!string_to_funcs.count(instruction.action)) {
            fprintf(stderr, "unrecognized action %s\n", instruction.action.c_str());
            exit(EXIT_FAILURE);
        }
        return (this->*(string_to_funcs.at(instruction.action).exec))(inSequences, instruction);
    }

};

#endif /* GFASTATS_SAK_H */
