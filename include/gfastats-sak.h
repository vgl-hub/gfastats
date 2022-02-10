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
    {"EXCISE", 3},
    {"REMOVE", 4},
    {"ERASE", 5},
    {"ADD", 6},
    {"REPLACE", 7},
    {"INVERT", 8},
    {"RVCP", 9}
};

struct Instruction {
    
    std::string action;
    std::string scaffold1;
    std::string scaffold2;
    std::string contig1;
    std::string contig2;
    std::string comment1;
    std::string comment2;
    std::string gHeader = "";
    
    char sId1Or, sId2Or;
    
    unsigned int dist, start, end;
    
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
                
                instruction.gHeader = arguments[4];
                
                instruction.scaffold1 = arguments[5];
                
                instruction.comment1 = arguments[6];
                
                break;
            }
                
            case 2: { // SPLIT
                
                instruction.contig1 = arguments[1];
                
                instruction.contig2 = arguments[2];
                
                break;
            }

            case 3: { // EXCISE
                
                instruction.contig1 = arguments[1];
                
                if (arguments[2] != "") {
                
                    instruction.dist = stoi(arguments[2]);
                    
                }else{
                    
                    instruction.dist = 0;
                    
                }
                
                instruction.gHeader = arguments[3];
                
                break;
            }
                
            case 4: { // REMOVE
                
                instruction.contig1 = arguments[1];
                
                break;
            }
                
            case 5: { // ERASE
                
                size_t pos1 = 0, pos2 = 0;
                
                pos1 = arguments[1].find(":");
                
                pos2 = arguments[1].find("-");
                
                instruction.contig1 = arguments[1].substr(0, pos1);
                
                instruction.start = stoi(arguments[1].substr(pos1+1, pos2));
                
                instruction.end = stoi(arguments[1].substr(pos2+1, arguments[1].size()+1));
                
                break;
            }
            
            default:
                fprintf(stderr, "unrecognized action %s\n", instruction.action.c_str());
                exit(1);
        }
        
        verbose(verbose_flag, "Instruction read");
        
        return instruction;
        
    }
    
    bool executeInstruction(InSequences& inSequences, Instruction instruction) {
        
        switch (string_to_case.count(instruction.action) ? string_to_case.at(instruction.action) : 0) {
            case 1: { // JOIN
                
                join(inSequences, instruction);
                
                break;
                
            }
                
            case 2: { // SPLIT
                
                split(inSequences, instruction);
                
                break;
                
            }
                
            case 3: { // EXCISE
                
                excise(inSequences, instruction);
                
                break;
                
            }

            case 4: { // REMOVE
                
                remove(inSequences, instruction);
                
                break;
                
            }
                
            case 5: { // ERASE
                
                erase(inSequences, instruction);
                
                break;
                
            }
                
            default:
                fprintf(stderr, "unrecognized action %s\n", instruction.action.c_str());
                return EXIT_FAILURE;
        }
        
        return true;
        
    }
    
    bool join(InSequences& inSequences, Instruction instruction) { // joins two sequences via a gap based on instruction
        
        InGap gap;
                
        gap.newGap(inSequences.gapUniqN+1, inSequences.headersToIds[instruction.contig1], inSequences.headersToIds[instruction.contig2], instruction.sId1Or, instruction.sId2Or, instruction.dist, instruction.gHeader); // define the new gap
            
        inSequences.gapUniqN++;
        
        inSequences.appendGap(gap); // introduce the new gap
        
        return true;
        
    }
    
    bool split(InSequences& inSequences, Instruction instruction) { // splits two sequences removing the gap in between based on instruction
        
        inSequences.removeGaps(&instruction.contig1, &instruction.contig2); // remove the gap
        
        return true;
        
    }

    bool excise(InSequences& inSequences, Instruction instruction) { // removes a sequence, removing also edges if present
        
        std::vector<InGap> oldGaps = inSequences.getGap(&instruction.contig1); // get neighbour gaps
        
        inSequences.removeGaps(&instruction.contig1); // remove the gaps associated with the excised contig
        
        if (instruction.dist > 0 && oldGaps[0].getsId1() != oldGaps[0].getsId2() && oldGaps[1].getsId1() != oldGaps[1].getsId2()) { // terminal gaps are not allowed to create new gaps when excised
        
            InGap gap;
        
            gap.newGap(inSequences.gapUniqN+1, oldGaps[0].getsId1(), oldGaps[1].getsId2(), oldGaps[0].getsId1Or(), oldGaps[1].getsId2Or(), instruction.dist, instruction.gHeader); // define the new gap
            
            inSequences.gapUniqN++;
            
            inSequences.appendGap(gap); // introduce the new gap
        
        }
        
        inSequences.updateGapLens();
            
        return true;
        
    }
    
    bool remove(InSequences& inSequences, Instruction instruction) { // removes a sequence, removing also edges if present
        
        inSequences.removeGaps(&instruction.contig1); // remove the gaps associated with contig1
        
        inSequences.removeSegment(&instruction.contig1); // remove the segment
        
        inSequences.updateGapLens();
        
        return true;
        
    }
    
    bool erase(InSequences& inSequences, Instruction instruction) { // removes a sequence, removing also edges if present
        
        inSequences.inSegments[inSequences.headersToIds[instruction.contig1]].trimSegment(instruction.start, instruction.end); // trim segment
        
        inSequences.changeTotSegmentLen(instruction.start-instruction.end);
        
        return true;
        
    }
    
};

#endif /* GFASTATS_SAK_H */
