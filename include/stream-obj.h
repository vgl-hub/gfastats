#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

class StreamObj {
    
    std::streambuf* buffer;
    std::shared_ptr<std::istream> stream;
    std::ifstream ifs;
    zstream::igzstream zfin, zin;
    bool file = false, gzip = false;
//    std::istringstream strm;
    char* content = new char[10000];
    std::condition_variable mutexCondition;
//    bool decompress = true, done = false;
    
public:
    
    StreamObj() :
    zfin(ifs), zin(std::cin) {}
    
    ~StreamObj(){this->closeStream(); delete[] content;}
    
    bool isGzip(std::streambuf* buffer);
    
//    void decompressBuf(std::streambuf* buffer);
//
//    void readBuf(std::streambuf* buffer);
    
    std::shared_ptr<std::istream> openStream(UserInput &userInput, char type);
    
    void closeStream();
    
    std::string type();
    
    std::shared_ptr<std::istream> returnStream();
    
};

#endif /* STREAM_OBJ_H */
