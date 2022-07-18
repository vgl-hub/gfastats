//
//  gfastats.h
//
//
//  Created by Giulio Formenti on 12/30/21.
//

#ifndef GFASTATS_H
#define GFASTATS_H

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <fstream>
#include <sstream>

#include <unistd.h>
#include <getopt.h>

#include <vector>  //required for zstream
#include <stack>
#include <queue>
#include <string.h>
#include <algorithm> //required for zstream
#include <cstring> //required for zstream
#include <tuple> // for graph manipulation
#include <cctype> // toupper()
#include <iomanip>

#include <chrono>
#include <memory>

#include <thread>
#include <mutex>
#include <condition_variable>

#include "gfastats-log.h"
Log lg;

#include "uid-generator.h"

#include "bed.h"

#include "gfastats-global.h" // global variables
#include "gfastats-struct.h"
#include "gfastats-functions.h" // global functions

#include "threadpool.h"

#include <parallel_hashmap/phmap.h>

#include "zlib.h"
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>
#include <zstream/ozstream.hpp>
#include <zstream/ozstream_impl.hpp>

#include "gfa-lines.h"

#include "gfa.h" // gfa classes
#include "gfastats-sak.h" // swiss army knife classes

#include "gfastats-output.h" // output classes
#include "input.h"

#endif /* GFASTATS_H */
