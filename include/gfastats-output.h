//
//  gfastats-classes.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef gfastatsoutput_h
#define gfastatsoutput_h

//classes
class Report {
private:
    unsigned int counter = 0;
    
public:
    bool generateSeqReport (InSequences &inSequences, InSequence &inSequence, int &outSequence_flag) {
        
        while (counter < inSequences.getScaffN()) {
            
            inSequence = inSequences.getInSequence(counter);
            
            std::cout<<output("Seq")<<counter+1<<std::endl;
            std::cout<<output("Header")<<inSequence.getSeqHeader()<<std::endl;
            std::cout<<output("Comment")<<inSequence.getSeqComment()<<std::endl;
            std::cout<<output("Total sequence length")<<inSequence.getSeqScaffLen()<<std::endl;
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
                std::cout<<output("Quality")<<inSequence.getInSequenceQuality()<<std::endl;
                
            }
            
            std::cout<<std::endl;
            counter++;
            
        }
        
        counter = 0;
        
        return true;
        
    }
    
    bool outFile (InSequences &inSequences, InSequence &inSequence, int splitLength, std::string &outSeq) {
        
        const static std::unordered_map<std::string,int> string_to_case{
            {"fasta",1},
            {"fa",1},
            {"fasta.gz",1},
            {"fa.gz",1},
            {"fastq",2},
            {"fq",2},
            {"fastq.gz",2},
            {"fq.gz",2},
            {"gfa",3},
            {"gfa.gz",3}
        };
        
        bool gzip = false;
        bool outFile = false;
        
        std::string path = rmFileExt(outSeq);
        std::string ext = getFileExt(outSeq);
        
        if(ext == "gz") {
            
            ext = getFileExt(path) + ".gz";
            path = rmFileExt(path);
            gzip = true;
            
        }
        
        if (string_to_case.find(path) == string_to_case.end()) {
            
            outFile = true;
            
        }else{
            
            ext = outSeq;
            
        }
        
        std::string output;
        
        std::unique_ptr<std::ostream> stream;
        
        std::ofstream os(outSeq);
        zstream::ogzstream zout(os);
        
        if (gzip && outFile) {
            
            stream = make_unique<std::ostream>(zout.rdbuf());
            
        }else if (!gzip && outFile){
            
            stream = make_unique<std::ostream>(os.rdbuf());
            
        }else{
            
            os.close();
            remove(outSeq.c_str());
            
            stream = make_unique<std::ostream>(std::cout.rdbuf());
            
        }
        
        switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
                
            case 1: {
                while (counter < inSequences.getScaffN()) {
                    
                    inSequence = inSequences.getInSequence(counter);
                    
                    *stream<<">"<<inSequence.getSeqHeader()<<" "<<inSequence.getSeqComment()<<std::endl;
                    
                    if (splitLength != 0) {
                        
                        textWrap(inSequence.getInSequence(), *stream, splitLength);
                        
                    }else{
                        
                        output = inSequence.getInSequence();
                        
                    }
                    
                    *stream<<output<<std::endl;
                    output = "";
                    
                    counter++;
                    
                }
                
            }
                
            case 2: {
                while (counter < inSequences.getScaffN()) {
                    
                    inSequence = inSequences.getInSequence(counter);
                    
                    *stream<<"@"<<inSequence.getSeqHeader()<<" "<<inSequence.getSeqComment()<<"\n"<<inSequence.getInSequence()<<"\n+\n"<<inSequence.getInSequenceQuality()<<"\n";
                    
                    counter++;
                    
                }
                
            }
                
            //case 0: {std::cout<<"Unrecognized output: "<<outSeq;}//undefined case
                
        }
        
        return true;
        
    }
    
};


#endif /* gfastats-output_h */
