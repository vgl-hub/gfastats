//
//  gfastats-global.h
//  
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFASTATS_GLOBAL_H
#define GFASTATS_GLOBAL_H

//global
static auto start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

// flags are global variables
static short int tabular_flag;
static int verbose_flag;
static int seqReport_flag;
static int outSequence_flag;
static int nstarReport_flag;
static int outSize_flag;
static int outCoord_flag;
static int outFile_flag;
static int stats_flag;
static int cmd_flag;
static int rmGaps_flag;
static int hc_compress_flag;

#endif /* GFASTATS_GLOBAL_H */
