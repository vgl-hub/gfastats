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
#include <string.h>
#include <algorithm> //required for zstream
#include <cstring> //required for zstream
#include <tuple> // for graph manipulation
#include <cctype> // toupper()
#include <iomanip>
#include <numeric>

#include <chrono>
#include <memory>

#include <parallel_hashmap/phmap.h>

#include <zlib/zlib.h>
#include <zstream/zstream_common.hpp>
#include <zstream/izstream.hpp>
#include <zstream/izstream_impl.hpp>
#include <zstream/ozstream.hpp>
#include <zstream/ozstream_impl.hpp>

#include "gfastats-global.h" // global variables
#include "gfastats-functions.h" // global functions

#include "gfastats-gfa.h" // gfa classes
#include "gfastats-sak.h" // swiss army knife classes

#include "gfastats-input.h" // input classes
#include "gfastats-output.h" // output classes

#endif /* GFASTATS_H */
