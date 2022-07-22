#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>

#include <iostream>
#include <fstream>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "functions.h" // global functions

#include "log.h"
#include "global.h"
#include "uid-generator.h"

#include "gfa-lines.h"

#include "threadpool.h"
#include "gfa.h"
#include "sak.h" // swiss army knife

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "stream-obj.h"

#include "input-agp.h"


void readAgp (InSequences& inSequences, UserInput& userInput) {

    inSequences.updateStats();
    
    StreamObj streamObj;

    std::shared_ptr<std::istream> stream;

    std::string pHeaderNew, pHeader1, pHeader2, gHeader, instruction, coord1, coord2, line;
    char pId1Or = '+', pId2Or;

    unsigned int pUId = 0, pUId1 = 0, pUId2 = 0, gUId = 0, dist = 0, seqLen, pathLen, start1 = 0, end1 = 0, start2 = 0, end2 = 0;
    phmap::flat_hash_map<std::string, unsigned int>* hash;
    phmap::flat_hash_map<std::string, unsigned int>::const_iterator got;

    std::vector<std::string> arguments; // line arguments
    std::vector<unsigned int> oldPaths; // vector of paths flagged to be removed only at the end

    std::vector<InPath> paths = inSequences.getInPaths();

    for (InPath path : paths) {
        
        oldPaths.push_back(path.getpUId());
        
    }

    stream = streamObj.openStream(userInput, 'a');

    std::queue<std::string> nextLines;

    while (true) {
        if(nextLines.size() > 0) {
            line = nextLines.front();
            nextLines.pop();
        } else if(!getline(*stream, line)) {
            break;
        }
        
        std::istringstream iss(line); // line to string
        
        arguments = readDelimited(line, "\t", "#"); // read the columns in the line
        
        if (arguments.size() == 0) {continue;}
        
        pHeaderNew = arguments[0]; // this is the current header
        
        if (arguments[4] == "W") { // this is an old path
            
            if (!discoverPaths_flag) {
            
                pHeader1 = arguments[5];
            
            }else{
                
                pHeader1 = arguments[5] + "_path";
            }
            
            pId1Or = arguments[8][0];
            
            hash = inSequences.getHash1();
            
            got = hash->find(pHeader1); // get the headers to uIds table
            
            if (got != hash->end()) { // this is not the first time we see this path
                
                pUId1 = got->second;
                
            }else{
                
                fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader1.c_str()); // sequence not found
                
                continue;
                
            }
            
            pathLen = inSequences.pathLen(pUId1);
            
            start1 = stoi(arguments[6]);
            end1 = stoi(arguments[7]);
            
            seqLen = end1 - start1 + 1;
            
            if(seqLen != pathLen) {

                fprintf(stderr, "Warning: sequence length (%u) differs from path length (%u). Subsetting (%s).\n", seqLen, pathLen, pHeader1.c_str());

            }else{
                
                start1 = 0;
                end1 = 0;
                
            }
            
            getline(*stream, line);
            nextLines.push(line);
            
            arguments = readDelimited(line, "\t", "#"); // read the next sequence
            
            if(pHeaderNew != arguments[0]) { // if this path does not need to be joined to anything that follows, we create a new path
                
                InPath path;
                
                pUId = inSequences.getuId();
                
                path.newPath(pUId, pHeaderNew);
                
                std::vector<PathComponent> pathComponents = inSequences.getInPath(pUId1).getComponents();
                
                path.append({std::begin(pathComponents), std::end(pathComponents)});
                
                inSequences.insertHash(pHeaderNew, pUId);
                
                inSequences.uId.next();
                
                inSequences.addPath(path);
                
                if(seqLen != pathLen) { // if it also needs to be trimmed
                    
                    inSequences.trimPathByUId(pUId, start1, end1);
                    
                }
                
                if(pId1Or == '-') {
                    
                    inSequences.revComPath(pUId);
                    
                }
                
            }
        
            
        }else if(arguments[4] == "N" || arguments[4] == "U"){

            hash = inSequences.getHash1();
            
            got = hash->find(pHeader1); // get the headers to uIds table (remove sequence orientation in the gap first)
            
            if (got == hash->end()) { // this is the first time we see this path
                
                fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader1.c_str()); // if the preceding sequence was not found we do not introduce a gap
                
                continue;
                
            }
            
            gUId = inSequences.getuId();
            
            if (arguments[6] == "scaffold") {
                
                hash = inSequences.getHash1();
                
                got = hash->find("gap"+std::to_string(gUId)); // get the headers to uIds table
                
                while (got != hash->end()) { // this is not the first time we see this path
                    
                    gUId++;
                    
                    got = hash->find("gap"+std::to_string(gUId)); // get the headers to uIds table
                    
                    inSequences.uId.next();
                }
            
                inSequences.uId.next();
                gHeader = "gap"+std::to_string(gUId);
            
            }else{
                
                gHeader = arguments[6];
                inSequences.uId.next();
                
            }
            
            inSequences.insertHash(gHeader, gUId);
            
            
            dist = stoi(arguments[5]);
            
            getline(*stream, line);
            
            arguments = readDelimited(line, "\t", "#"); // read the next sequence
            
            if (arguments.size() == 0) {continue;}
            
            if (!discoverPaths_flag) {
            
                pHeader2 = arguments[5];
            
            }else{
                
                pHeader2 = arguments[5] + "_path";
            }
                
            pId2Or = arguments[8][0];
            
            hash = inSequences.getHash1();
            
            got = hash->find(pHeader2); // get the headers to uIds table (remove sequence orientation in the gap first)
            
            if (got != hash->end()) { // this is not the first time we see this path
                
                pUId2 = got->second;
                
            }else{
                
                fprintf(stderr, "Warning: sequence missing from the path set (%s). Skipping.\n", pHeader2.c_str()); // sequence not found
                
                continue;
                
            }
            
            pathLen = inSequences.pathLen(pUId2);
            
            start2 = stoi(arguments[6]);
            end2 = stoi(arguments[7]);
            
            seqLen = end2 - start2 + 1;
            
            if(seqLen != pathLen) {

                fprintf(stderr, "Warning: sequence length (%u) differs from path length (%u). Subsetting (%s).\n", seqLen, pathLen, pHeader2.c_str());

            }else{
                
                start2 = 0;
                end2 = 0;
                
            }
            
            SAK sak;
            
            coord1 = start1 != 0 ? "(" + std::to_string(start1) + ":" + std::to_string(end1) + ")" : "";
            coord2 = start2 != 0 ? "(" + std::to_string(start2) + ":" + std::to_string(end2) + ")" : "";
            
            instruction = "JOIN\t" + pHeader1 + coord1 + pId1Or + "\t" + pHeader2 + coord2 + pId2Or + "\t" + std::to_string(dist) + "\t" + gHeader + "\t" + pHeaderNew + "\t" + std::to_string(gUId);
            
            fprintf(stderr, "%s\n", instruction.c_str());
            
            sak.executeInstruction(inSequences, sak.readInstruction(instruction));
            
            pHeader1 = pHeaderNew;
            start1 = 0;
            end1 = 0;
            pId1Or = '+';

        }
        
    }
        
    for (unsigned int pUId : oldPaths) { // remove paths left
        
        inSequences.removePath(pUId, false, true); // silently remove the original paths that were not joined or duplicated
        
    }

}
