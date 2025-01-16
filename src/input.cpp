#include <stdlib.h>
#include <string>

#include <istream>
#include <fstream>
#include <sstream>

#include <parallel-hashmap/phmap.h>

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"

#include "gfa-lines.h"
#include "gfa.h"
#include "sak.h"

#include "stream-obj.h"

#include "input-agp.h"
#include "input-filters.h"
#include "input-gfa.h"
#include "input.h"

void Input::load(UserInputGfastats userInput) {
    
    this->userInput = userInput;
    
}
    
void Input::read(InSequences& inSequences) {
    
    if (userInput.inSequence.empty()) {return;}
    
    threadPool.init(maxThreads); // initialize threadpool

    if (!userInput.inSak.empty() || userInput.pipeType == 'k') {
        
        StreamObj streamObj;
        
        stream = streamObj.openStream(userInput, 'k');
        
        SAK sak; // create a new swiss army knife
        
        while (getline(*stream, line)) {
            
            std::istringstream iss(line);
            
            instructions.push_back(sak.readInstruction(line)); // use the swiss army knife to read the instruction
            
        }
        
        lg.verbose("Finished reading SAK instructions");
        
    }
    
    if (!userInput.inBedInclude.empty() || userInput.pipeType == 'i') {
        
        StreamObj streamObj;
        stream = streamObj.openStream(userInput, 'i');
        
        while (getline(*stream, line)) {
            
            uint64_t begin = 0, end = 0;
            std::istringstream iss(line);
            iss >> bedHeader >> begin >> end;
            userInput.bedIncludeList.pushCoordinates(bedHeader, begin, end);
        }
        lg.verbose("Finished reading BED include list");
    }
    
    BedCoordinates bedExcludeList;
    
    if (!userInput.inBedExclude.empty() || userInput.pipeType == 'e') {
        
        StreamObj streamObj;
        stream = streamObj.openStream(userInput, 'e');
        
        while (getline(*stream, line)) {
            
            uint64_t begin = 0, end = 0;
            std::istringstream iss(line);
            iss >> bedHeader >> begin >> end;
            
            bedExcludeList.pushCoordinates(bedHeader, begin, end);
        }
        lg.verbose("Finished reading BED exclude list");
    }
    
    if (!userInput.inSequence.empty() || userInput.pipeType == 'f') {
        
        StreamObj streamObj;
        
        stream = streamObj.openStream(userInput, 'f');
        
        if (stream) {
            
            switch (stream->peek()) {
                    
                case '>': {
                    
                    stream->get();
                    
                    while (getline(*stream, newLine)) {
                        
                        if(userInput.bedIncludeList.size() - bedExcludeList.size() != 0 && userInput.bedIncludeList.size() - bedExcludeList.size() == inSequences.getPathN()) { // we have all the sequences needed
                            lg.verbose("Found all sequences, stop streaming input");
                            break;
                        }
                        size_t spacePos = newLine.find(" ");
                        seqHeader = newLine.substr(0, spacePos);
                        if (spacePos != std::string::npos)
                            seqComment = newLine.substr(spacePos + 1);
                        
                        std::string* inSequence = new std::string;
                        
                        getline(*stream, *inSequence, '>');
                        
                        lg.verbose("Individual fasta sequence read");
                        
                        Sequence* sequence = includeExcludeSeq(seqHeader, seqComment, inSequence, userInput.bedIncludeList, bedExcludeList);
                        
                        if (sequence != NULL) {
                            
                            sequence->seqPos = seqPos; // remember the order
                            
                            inSequences.appendSequence(sequence, userInput.hc_cutoff);
                            
                            seqPos++;
                            
                        }
                        
                    }
                    
                    break;
                }
                case '@': {
                    
                    while (getline(*stream, newLine)) { // file input
                        
                        if(userInput.bedIncludeList.size() - bedExcludeList.size() != 0 && userInput.bedIncludeList.size() - bedExcludeList.size() == inSequences.getPathN()) { // we have all the sequences needed
                            lg.verbose("Found all sequences, stop streaming input");
                            break;
                        
                        }
                        newLine.erase(0, 1);
                        size_t spacePos = newLine.find(" ");
                        seqHeader = newLine.substr(0, spacePos);
                        if (spacePos != std::string::npos)
                            seqComment = newLine.substr(spacePos + 1);
                        
                        std::string* inSequence = new std::string;
                        getline(*stream, *inSequence);
                        
                        getline(*stream, newLine);
                        
                        std::string* inSequenceQuality = new std::string;
                        getline(*stream, *inSequenceQuality);

                        Sequence* sequence = includeExcludeSeq(seqHeader, seqComment, inSequence, userInput.bedIncludeList, bedExcludeList, inSequenceQuality);
                        
                        if (sequence != NULL) {
                            
                            sequence->seqPos = seqPos; // remember the order
                        
                            inSequences.appendSequence(sequence, userInput.hc_cutoff);
                            
                            seqPos++;
                            
                        }
                        
                    }
                    
                    break;
                    
                }
                default: {
                    
                    readGFA(inSequences, userInput, stream, &bedExcludeList);
                    
                }
                
            }
            
            lg.verbose("End of file");
                
        }else{

            fprintf(stderr, "Stream not successful: %s", userInput.inSequence.c_str());
            exit(1);

        }
        
    }
    
    jobWait(threadPool);
    
    inSequences.sortSegmentsByOriginal();
    
    if (userInput.rmGaps_flag)
        inSequences.removeTerminalGaps();
    
    if (userInput.extractContigs_flag) {
        
        inSequences.clearGaps();
        inSequences.clearPaths();
        
    }
    
    if (userInput.extractContigs_flag || userInput.discoverPaths_flag)
        inSequences.discoverPaths();
    
    if (userInput.terminalOvlLen != 0)
        inSequences.discoverTerminalOverlaps(userInput.terminalOvlLen);
    
    if (!instructions.empty()) {
        
        lg.verbose("\nStarted instruction execution");
    
        SAK sak; // create a new swiss army knife
        
        for (Instruction instruction : instructions) { // execute swiss army knife instructions
            
            sak.executeInstruction(inSequences, instruction);
            lg.verbose(instruction.action + " instruction executed");
            
        }
    
    }
    
    if (!userInput.inAgp.empty() || userInput.pipeType == 'a')
        readAgp(inSequences, userInput);
    
    if (userInput.sortType == "ascending") {
        inSequences.sortPathsByNameAscending();
    }else if (userInput.sortType == "descending") {
        inSequences.sortPathsByNameDescending();
    }else if (userInput.sortType == "largest") {
        inSequences.sortPathsBySize(0);
    }else if (userInput.sortType == "smallest") {
        inSequences.sortPathsBySize(1);
    }else if (userInput.sortType != "none" && ifFileExists(userInput.sortType.c_str())){
            
        stream = std::make_unique<std::ifstream>(std::ifstream(userInput.sortType));
        
        std::string header;
        std::vector<std::string> headerList;
        
        while (getline(*stream, line)) { // read the file to vector
            
            std::istringstream iss(line);
            iss >> header;
            
            headerList.push_back(header);
            
        }
        
        inSequences.sortPathsByList(headerList);
        
    }else if(userInput.inAgp.empty() && !(userInput.pipeType == 'a')){
        inSequences.sortPathsByOriginal();
    }
    
    inSequences.updateStats();
    
    threadPool.join();
    
}
