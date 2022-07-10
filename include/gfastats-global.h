//
//  gfastats-global.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFASTATS_GLOBAL_H
#define GFASTATS_GLOBAL_H

#include <stdlib.h>
#include <chrono>

//global time
extern std::chrono::high_resolution_clock::time_point start;

// flags are global variables
extern short int tabular_flag;
extern int verbose_flag;
extern int seqReport_flag;
extern int outSequence_flag;
extern int nstarReport_flag;
extern int outSize_flag;
extern int outCoord_flag;
extern int outFile_flag;
extern int outBubbles_flag;
extern int stats_flag;
extern int cmd_flag;
extern int rmGaps_flag;
extern int discoverPaths_flag;
extern int hc_flag;
extern int hc_cutoff;
extern int maxThreads;

extern Log lg;
extern std::mutex mtx;

#endif /* GFASTATS_GLOBAL_H */
