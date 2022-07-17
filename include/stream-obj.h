//
//  gfastats-functions.h
//
//
//  Created by Giulio Formenti on 07/16/22.
//

#ifndef STREAM_OBJ_H
#define STREAM_OBJ_H

class StreamObj {
    
    std::streambuf* buffer;
    std::shared_ptr<std::istream> stream;
    std::ifstream ifs;
    zstream::igzstream zfin, zin;
    bool file = false, gzip = false;

public:
    
    StreamObj() :
    zfin(ifs), zin(std::cin) {}
    
    ~StreamObj(){this->closeStream();}
    
    bool isGzip(std::streambuf* buffer);
    std::shared_ptr<std::istream> openStream(UserInput &userInput, char type);
    void closeStream();
    std::string type();
    std::shared_ptr<std::istream> returnStream();
    
};

#endif /* STREAM_OBJ_H */
