//
//  gfastats-sak.h
//
//
//  Created by Giulio Formenti on 1/30/22.
//

#ifndef GFASTATS_SAK_H
#define GFASTATS_SAK_H

// unordered map to handle out correspondence in following switch statement
const static std::unordered_map<std::string,int> string_to_case{
    {"JOIN", 1},
    {"SPLIT", 2},
    {"REMOVE", 3},
    {"DELETE", 4},
    {"ADD", 5},
    {"REPLACE", 6},
    {"EXCISE", 7},
    {"INVERT", 8},
    {"RVCP", 9},
};

struct Instruction {
    
    std::string action;
    std::string scaffold1;
    std::string scaffold2;
    std::string contig1;
    std::string contig2;
    std::string comment1;
    std::string comment2;
    
    char sId1Or, sId2Or;
    
    unsigned int dist;
    
};

class SAK { // the swiss army knife
private:
    InSequences inSequences;
    InSegment inSegment1, inSegment2, inSegmentNew;
    std::string sId1Header, sId2Header;
    
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
        
        switch (string_to_case.count(instruction.action) ? string_to_case.at(instruction.action) : 0) {
            
            case 1: { // JOIN
                
                instruction.sId1Or = arguments[1].back(); // get sequence orientation in the gap
                
                arguments[1].pop_back(); // remove sequence orientation in the gap
                
                instruction.contig1 = arguments[1];
                
                instruction.sId2Or = arguments[2].back(); // get sequence orientation in the gap
                
                arguments[2].pop_back(); // remove sequence orientation in the gap
                
                instruction.contig2 = arguments[2];
                
                instruction.dist = stoi(arguments[3]);
                
                instruction.scaffold1 = arguments[4];
                
                instruction.comment1 = arguments[5];
                
                break;
            }
            
            default:
                fprintf(stderr, "unrecognized action %s\n", instruction.action.c_str());
                exit(1);
        }
        
        return instruction;
        
    }
    
    bool executeInstruction(InSequences& inSequences, Instruction instruction) {
        
        switch (string_to_case.count(instruction.action) ? string_to_case.at(instruction.action) : 0) {
            case 1: { // JOIN
                
                joinByGap(inSequences, instruction);
                
                break;
                
            }
            
            default:
                fprintf(stderr, "unrecognized action %s\n", instruction.action.c_str());
                return EXIT_FAILURE;
        }
        
        return true;
        
    }
    
    bool joinByGap(InSequences& inSequences, Instruction instruction) { // joins two sequences via a gap based on instruction in gfa format
        
        InGap gap;
        
        gap.newGap(inSequences.gapUniqN+1, inSequences.headersToIds[instruction.contig1], inSequences.headersToIds[instruction.contig2], instruction.sId1Or, instruction.sId2Or, instruction.dist); // define the new gap
        
        
//        if (instruction.scaffold1 != "") {
//
//            inSequences.inSegments[inSequences.headersToIds[instruction.contig1]].seqHeader = instruction.scaffold1;
//
//        }
//
//        if (instruction.comment1 != "") {
//
//            inSequences.inSegments[inSequences.headersToIds[instruction.contig1]].seqComment = instruction.comment1;
//
//        }
        
        inSequences.gapUniqN++;
        
        inSequences.appendGap(gap); // introduce the new gap
        
        return true;
        
    }
    
};

#endif /* GFASTATS_SAK_H */
