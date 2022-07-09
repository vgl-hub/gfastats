#ifndef GFASTATS_LOG_H
#define GFASTATS_LOG_H

struct Log {

    std::string log;
    unsigned int jobId;
    
    void verbose(std::string msg, bool overwrite = false); // verbose decorated output
    
    void add(std::string msg); // verbose decorated output
    
    void print();
    
    void setId (unsigned int i);
    
};

#endif /* GFASTATS_LOG_H */
