#include <memory>
#include <chrono>
#include <cstring>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <vector>
#include <iostream>
#include <unistd.h>
#include <string>

#include <gfastats-log.h>
#include <gfastats-global.h>
#include "gfastats-struct.h"

double elapsedTime(){ // compute runtime in verbose mode
    
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    start = std::chrono::high_resolution_clock::now();
    
    return elapsed.count();
    
}

//functions
bool checkTag(const char tag1[2], std::string tag2) {
    return tag1 == tag2;
}

bool isInt(const std::string &str) {
    return !str.empty() && str.find_first_not_of("0123456789") == std::string::npos;
}

double gfa_round(double d, uint32_t to=2) {
    if(std::isnan(d)) return NAN;
    unsigned int n=1;
    for(; to>0; to--) n *= 10;
    return std::round(d*n)/n;
}

std::vector<unsigned int> intervalSizes(std::vector<unsigned int> &intervalVec){ // compute sizes of a vector of intervals to be read as paired coordinates in 0-bed format
    
    std::vector<unsigned int> intervalVecLens;
    intervalVecLens.reserve(200);
    
    if (!intervalVec.empty()){
        std::vector<unsigned int>::const_iterator end = intervalVec.cend();
        
        for (std::vector<unsigned int>::const_iterator it = intervalVec.cbegin(); it != end;) {
            
            intervalVecLens.push_back(*(it+1) - *it); // compute size of the interval
            
            it = it + 2; // jump to next pair
            
        }
        
    }
    
    return intervalVecLens;
    
}

std::string output(std::string output){ // use tab delimiter if tabular flag is true
    
    if (tabular_flag) {
        
        output = output + "\t";
        
    }else{
        
        output = output + ": ";
        
    }
    
    return output;
}

bool isDash(char * optarg) { // check if user input is dash (substitute of input from pipe)
    
    return (strcmp(optarg, "-") == 0) ? true : false;
    
}

bool ifFileExists(const char * optarg) { // check if file exists
    
    if (!access (optarg, F_OK)) {
        
        return optarg;
        
    }else{
        
        std::cout<<"Error - file does not exist: "<<optarg<<std::endl;
        exit(1);
        
    }
    
}

void textWrap(std::string input, std::ostream& output, int width) { // generic text wrapper (useful for fasta output)
    
    std::string tmp;
    char cur = '\0';
    int i = 0;
    
    std::stringstream ss(input);
    
    while (ss.get(cur)) {
        if (++i == width+1) {

            output << tmp << '\n';
            i = tmp.length();
            tmp.clear();
            
        }else{
            
            output << tmp;
            tmp.clear();
            
        }
        
        tmp += cur;
        
    }
    
    output << tmp << '\n';
    tmp.clear();
    
}

std::string rmFileExt(const std::string& path) { // utility to strip file extension from file
    if (path == "." || path == "..")
        return path;

    size_t pos = path.find_last_of("\\/.");
    if (pos != std::string::npos && path[pos] == '.')
        return path.substr(0, pos);

    return path;
}

std::string getFileExt(const std::string& FileName) // utility to get file extension
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}

std::string revCom(std::string seq) { // reverse complement
    auto lambda = [](const char c) {
        switch (c) {
        case '*':
            return '*';
        case 'A':
            return 'T';
        case 'G':
            return 'C';
        case 'C':
            return 'G';
        case 'T':
            return 'A';
        case 'a':
            return 't';
        case 'g':
            return 'c';
        case 'c':
            return 'g';
        case 't':
            return 'a';
        case 'N':
        case 'n':
        case 'X':
        case 'x':
            return c;
        default:
            throw std::domain_error("Invalid nucleotide.");
        }
    };

    std::transform(seq.cbegin(), seq.cend(), seq.begin(), lambda);
    reverse(seq.begin(), seq.end());
    return seq;
}

std::string rev(std::string seq) { // reverse string
    
    reverse(seq.begin(), seq.end());
 
    return seq;
    
}

std::vector<std::string> readDelimited(std::string line, std::string delimiter, std::string skipLine = "") { // read line delimited by specific character, optionally skip lines starting with specific string

    std::vector<std::string> arguments;

    size_t pos = 0;
    
    if (skipLine != "" && line.substr(0, skipLine.size()) == skipLine) {
        
        return arguments;
        
    }

    while ((pos = line.find(delimiter)) != std::string::npos) {
        
        arguments.push_back(line.substr(0, pos));
        
        line.erase(0, pos + delimiter.length());
            
    }
    
    arguments.push_back(line); // last column
        
    return arguments;
    
}

bool isNumber(const std::string& str)
{
    for (char const &c : str) {
        if (std::isdigit(c) == 0) return false;
    }
    return true;
}

void revComPathComponents(std::vector<PathComponent>& pathComponents) {
    
    for (PathComponent& component : pathComponents) {
        
        if (component.orientation != '0') {
        
            component.orientation = (component.orientation == '+' ? '-' : '+');
        
        }
        
    }
    
    std::reverse(pathComponents.begin(), pathComponents.end());
    
}

// bed coords are bed coords of compressed sequence
void homopolymerCompress(std::string *sequence, std::vector<std::pair<unsigned long long int, unsigned long long int>> &bedCoords, unsigned int cutoff) {
    unsigned int index=0, length, new_length=0;

    auto lambda = [&length, &index, &bedCoords, &sequence, &new_length, &cutoff](int i){
        length = i-index;
        if(length > cutoff) {
            bedCoords.push_back({new_length, new_length+length});
        }
        int num = length > cutoff ? 1 : length;
        memset(&((*sequence)[new_length]), (*sequence)[index], num);
        new_length += num;
    };

    for(unsigned int i=1; i<sequence->length(); ++i) {
        if((*sequence)[i] == (*sequence)[index]) continue;
        lambda(i);
        index = i;
    }
    lambda(sequence->length());
    sequence->resize(new_length);
}

// bed coords are bed coords of compressed sequence
void homopolymerDecompress(std::string *sequence, const std::vector<std::pair<unsigned long long int, unsigned long long int>> &bedCoords) {
    std::string ret="";
    ret.reserve(sequence->length()*2); // random guess for final sequence length to reduce resizes
    for(unsigned int i=0, ci=0, len; i<sequence->length(); ++i) {
        if(ci < bedCoords.size() && i == bedCoords[ci].first) {
            len = bedCoords[ci].second - bedCoords[ci].first;
            ++ci;
        } else {
            len = 1;
        }
        ret += std::string(len, (*sequence)[i]);
    }
    ret.shrink_to_fit();
    *sequence = ret;
}

unsigned int homopolymerRunsCount(const std::string &sequence, unsigned int threshhold) {
    unsigned int runs = 0;
    unsigned int currentRun=0;
    char prev=0;
    for(const char &curr : sequence) {
        if(prev != curr) {
            if(currentRun >= threshhold) {
                ++runs;
            }
            currentRun = 0;
        }
        else {
            ++currentRun;
        }
        prev = curr;
    }
    if(currentRun >= threshhold) { // loop wont catch run at end of sequence
        ++runs;
    }
    return runs;
}

// bed coords of uncompressed sequence
void homopolymerBedCoords(std::string *sequence, std::vector<std::pair<unsigned int, unsigned int>> &bedCoords, unsigned int cutoff) {
    int index = 0;
    for(unsigned int i=1; i < sequence->size(); ++i) {
        if((*sequence)[i] == (*sequence)[i-1]) continue;
        if(i-index > cutoff) {
            bedCoords.push_back({index, i});
        }
        index = i;
    }
    if((*sequence)[sequence->size()-1] == (*sequence)[sequence->size()-2] && sequence->size()-index > cutoff) {
        bedCoords.push_back({index, sequence->size()});
    }
}

void traverseInSequence(Sequence sequence);
