//
//  gfastats-functions.h
//
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFASTATS_FUNCTIONS_H
#define GFASTATS_FUNCTIONS_H

//typedef
typedef std::tuple<char, unsigned int, char, unsigned int, unsigned int> Tuple; // tuple for gap edges orientation|segment_id|orientation|dist|edge_id
typedef std::tuple<char, unsigned int, char> EdgeTuple; // tuple for gap edges orientation|id|orientation
typedef std::tuple<char, unsigned int, char> PathTuple; // tuple for paths type|id|orientation

//templates
template<typename T, typename... Args> // unique pointer to handle different types of istreams and ostreams
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

//functions
bool isInt(const std::string &str) {
    return !str.empty() && str.find_first_not_of("0123456789") == std::string::npos;
}

double elapsedTime(){ // compute runtime in verbose mode
    
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    start = std::chrono::high_resolution_clock::now();
    
    return elapsed.count();
    
}

void verbose(std::string msg) { // verbose decorated output
    
    if(verbose_flag) {
        
        std::cout << msg<< " (done in " << elapsedTime() << " s).\n";
        
        elapsedTime();
        
    }
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

bool determineGzip(std::string iFastaFileArg) { // check first two bytes for gzip compression
    
    std::ifstream stream(iFastaFileArg);
    
    unsigned char buffer[2];
    stream.read((char*)(&buffer[0]), 2);
    
    stream.clear();
    stream.seekg(0, stream.beg);
    
    if (buffer[0] == 0x1f && (buffer[1] == 0x8b)) {
        
        stream.close();
        
        return true;
        
    }else{
        
        return false;
        
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

#endif /* GFASTATS_FUNCTIONS_H */
