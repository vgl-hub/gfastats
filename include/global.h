#ifndef GLOBAL_H
#define GLOBAL_H

#include <mutex>
#include <chrono>
#include <queue>
#include <thread>
#include <functional>

#include "log.h"
#include "threadpool.h"

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
extern int extractContigs_flag;
extern int hc_flag;
extern int hc_cutoff;
extern int terminalOvlLen;
extern int maxThreads;

extern Log lg;
extern std::mutex mtx;
extern ThreadPool<std::function<bool()>> threadPool;

#endif /* GLOBAL_H */
