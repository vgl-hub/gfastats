//
//  gfastats-classes.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef gfastatsoutputClasses_h
#define gfastatsoutputClasses_h

//classes
class Report {
private:
    unsigned int counter = 0;
    
public:
    bool generateSeqReport (InSequences &inSequences, InSequence &inSequence, int &outSequence_flag) {
        
        while (counter < inSequences.getScaffN()) {
            
            inSequence = inSequences.getInSequence(counter);
            
            std::cout<<output("Seq")<<counter+1<<std::endl;
            std::cout<<output("Header")<<inSequence.getFastaHeader()<<std::endl;
            std::cout<<output("Comment")<<inSequence.getFastaComment()<<std::endl;
            std::cout<<output("Total sequence length")<<inSequence.getFastaScaffLen()<<std::endl;
            std::cout<<output("Total contig length")<<inSequence.getContigSum()<<std::endl;
            std::cout<<output("# contig")<<inSequence.getContigN()<<std::endl;
            std::cout<<output("Total gap length")<<inSequence.getGapSum()<<std::endl;
            std::cout<<output("# gaps")<<inSequence.getGapN()<<std::endl;
            
            printf("%s%u, %u, %u, %u\n",output("Base composition (ACGT)").c_str(), inSequence.getA(),
                   inSequence.getC(),
                   inSequence.getG(),
                   inSequence.getT());
            printf("%s%.2f\n",output("GC content %").c_str(), inSequence.computeGCcontent());
            std::cout<<output("# soft-masked bases")<<inSequence.getLowerCount()<<std::endl;
            
            
            if (outSequence_flag) {
                
                std::cout<<output("Sequence")<<inSequence.getInSequence()<<std::endl;
                
            }
            
            std::cout<<std::endl;
            counter++;
            
        }
        
        counter = 0;
        
        return true;
        
    }
    
    bool outFasta (InSequences &inSequences, InSequence &inSequence, short unsigned int splitLength, int &outSequence_flag) {
        
        while (counter < inSequences.getScaffN()) {
            
            inSequence = inSequences.getInSequence(counter);
            
            std::cout<<">"<<inSequence.getFastaHeader()<<" "<<inSequence.getFastaComment()<<std::endl;
            
            if (splitLength != 0) {
                
                unsigned int pos = 0;
                std::string line;
                
                for (char& base : inSequence.getInSequence())
                {
                    
                    line += base;
                    
                    if (pos == splitLength) {
                        
                        std::cout<<line;
                        std::cout<<std::endl;
                        
                        line = "";
                        pos = 0;
                        
                    }
                    
                    pos++;
                    
                }
                
                if (inSequence.getInSequence().length() % splitLength != 0) {
                    
                    std::cout<<std::endl;
                    
                }
                
            }else{
                
                std::cout<<inSequence.getInSequence()<<std::endl;
                
            }
            
            counter++;
            
        }
        
        return true;
        
    }
    
};


#endif /* gfastats-classes_h */
