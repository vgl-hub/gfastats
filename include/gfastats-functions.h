//
//  gfastats-functions.h
//
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef gfastatsFunctions_h
#define gfastatsFunctions_h

//functions

bool isInt(const std::string &str) {
  return !str.empty() && str.find_first_not_of("0123456789") == std::string::npos;
}

double elapsedTime(){
    
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    start = std::chrono::high_resolution_clock::now();
    
    return elapsed.count();
    
}

void verbose(int &verbose_flag, std::string msg) {
    
    if(verbose_flag) {
        
        std::cout << msg<< " (done in " << elapsedTime() << " s).\n";
        
        elapsedTime();
        
    }
}

std::vector<unsigned int> bedIntervalSizes(std::vector<unsigned int> &intervalVec){
    
    std::vector<unsigned int> intervalVecLens;
    intervalVecLens.reserve(200);
    
    if (!intervalVec.empty()){
        std::vector<unsigned int>::const_iterator end = intervalVec.cend();
        
        for (std::vector<unsigned int>::const_iterator it = intervalVec.cbegin(); it != end;) {
            
            intervalVecLens.push_back(*(it+1) - *it);
            
            it = it + 2;
            
        }
        
    }
    
    return intervalVecLens;
    
}

std::string output(std::string output){
    
    if (tabular_flag) {
        
        output = output + "\t";
        
    }else{
        
        output = output + ": ";
        
    }
    
    return output;
}

bool loadBinaryFile(const std::string& filename, std::string& contents);
bool gzipInflate(const std::string& compressedBytes, std::string& uncompressedBytes);

std::string loadGzip(std::string iFastaFileArg) {

    std::string fileData;
    if (!loadBinaryFile(iFastaFileArg, fileData)) {
        printf("Error loading input file: %s", iFastaFileArg.c_str());
    }

    std::string data;
    if (!gzipInflate(fileData, data)) {
        printf("Error decompressing input file: %s", iFastaFileArg.c_str());
        
    }
    
    return(data);
    
}

/* Reads a file into memory. */
bool loadBinaryFile(const std::string& filename, std::string& contents) {
    // Open the gzip file in binary mode
    FILE* f = fopen(filename.c_str(), "rb");
    if (f == NULL)
        return false ;
    
    // Clear existing bytes in output std::vector
    contents.clear();
    
    // Read all the bytes in the file
    int c = fgetc(f);
    while (c != EOF) {
        contents +=  (char) c ;
        c = fgetc(f);
    }
    fclose (f);
    
    return true ;
}

//zlib
bool gzipInflate(const std::string& compressedBytes, std::string& uncompressedBytes) {
    if (compressedBytes.size() == 0) {
        uncompressedBytes = compressedBytes ;
        return true ;
    }
    
    uncompressedBytes.clear() ;
    
    unsigned full_length = compressedBytes.size() ;
    unsigned half_length = compressedBytes.size() / 2;
    
    unsigned uncompLength = full_length ;
    char* uncomp = (char*) calloc(sizeof(char), uncompLength);
    
    z_stream strm;
    strm.next_in = (Bytef *) compressedBytes.c_str();
    strm.avail_in = compressedBytes.size() ;
    strm.total_out = 0;
    strm.zalloc = Z_NULL;
    strm.zfree = Z_NULL;
    
    bool done = false ;
    
    if (inflateInit2(&strm, (16+MAX_WBITS)) != Z_OK) {
        free(uncomp);
        return false;
    }
    
    while (!done) {
        // If our output buffer is too small
        if (strm.total_out >= uncompLength) {
            // Increase size of output buffer
            char* uncomp2 = (char*) calloc(sizeof(char), uncompLength + half_length);
            memcpy(uncomp2, uncomp, uncompLength);
            uncompLength += half_length ;
            free(uncomp);
            uncomp = uncomp2 ;
        }
        
        strm.next_out = (Bytef *) (uncomp + strm.total_out);
        strm.avail_out = uncompLength - strm.total_out;
        
        // Inflate another chunk.
        int err = inflate (&strm, Z_SYNC_FLUSH);
        if (err == Z_STREAM_END) done = true;
        else if (err != Z_OK)  {
            break;
        }
    }
    
    if (inflateEnd (&strm) != Z_OK) {
        free(uncomp);
        return false;
    }
    
    for (size_t i=0; i<strm.total_out; ++i) {
        uncompressedBytes += uncomp[i];
    }
    free(uncomp);
    return true ;
}

bool isDash(char * optarg) {
    
    return (strcmp(optarg, "-") == 0) ? true : false;

}

bool ifFileExists(char * optarg) {

    if (!access (optarg, F_OK)) {
        
        return optarg;
        
    }else{
        
        std::cout<<"Error - file does not exist: "<<optarg<<std::endl;
        exit(1);
        
    }

}

bool determineGzip(std::string iFastaFileArg) {
    
    std::ifstream stream(iFastaFileArg);
    
    unsigned char buffer[2];
    stream.read((char*)(&buffer[0]), 2) ;
    
    stream.clear();
    stream.seekg(0, stream.beg);
    
    if (buffer[0] == 0x1f && (buffer[1] == 0x8b)) {
        
        stream.close();
        
        return true;
        
    }else{
        
        return false;
            
    }
    
}

//templates
template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

#endif /* gfastats-Functions_h */
