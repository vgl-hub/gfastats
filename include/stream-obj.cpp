#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <vector>

#include "bed.h"
#include "gfastats-struct.h"

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>

#include "gfastats-log.h"
extern Log lg;

#include "stream-obj.h"



std::string StreamObj::type() {
    
    std::string type;
    type =  file ? "file" : "pipe";
    type += "/";
    type += gzip ? "gzip" : "plain";
    
    return type;
    
}

bool StreamObj::isGzip(std::streambuf* buffer) {
    
    unsigned char header[2];
    
    header[0] = buffer->sgetc();
    header[1] = buffer->snextc();
    buffer->sungetc();
    
    return (header[0] == 0x1f && header[1] == 0x8b) ? true : false;
    
}

std::shared_ptr<std::istream> StreamObj::openStream(UserInput &userInput, char type) {
    
    file = userInput.pipeType != type ? true : false;
    
    if (file) { // input is from file
        
        ifs.close();
        
        ifs.open(userInput.file(type)); // this stream takes input from a plain file

        buffer = ifs.rdbuf();
        
        gzip = isGzip(buffer);

        if (gzip) {
            
            zfin.open();

            buffer = zfin.rdbuf();

        }

    }else{
        
        buffer = std::cin.rdbuf();
        
        gzip = isGzip(buffer);

        if (gzip) {
            
            zin.open();
            
            buffer = zin.rdbuf();

        }
        
    }
    
    return std::make_shared<std::istream>(buffer);
    
}

std::shared_ptr<std::istream> StreamObj::returnStream() {
    
    return stream;
    
}

void StreamObj::closeStream() {

    if (gzip) {

        zfin.read_footer();

        if (zfin.check_crc()) {

            lg.verbose("Crc check successful");

        }else{

            lg.verbose("Warning: crc check unsuccessful. Check input file");

        }

    }

}
