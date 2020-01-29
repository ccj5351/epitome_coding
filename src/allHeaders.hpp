#ifndef _ALLMYHEADERS_HPP
#define _ALLMYHEADERS_HPP

// see http://stackoverflow.com/questions/14448433/error-c4996-fopen-not-declared
#define _CRT_SECURE_NO_DEPRECATE
#include <cstdio> // C header

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <complex> // C++, not C
#include <fftw3.h>
#include <math.h>
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <cfloat> // DBL_MAX = 1E+37 or greater
#include <cstdint> // This header defines a set of integral type aliases with specific width requirements, 
                   // along with macros specifying their limits and macro functions to create values of these types.

#include <random>  // a random number generator
#include <tuple> // e.g., std::tuple, std::get, std::tie, std::ignore

// boost filesystem
#define BOOST_FILESYSTEM_NO_DEPRECATED // For new code, defining "BOOST_FILESYSTEM_NO_DEPRECATED" before including filesystem headers is strongly recommended. 
                                       // This prevents inadvertent use of old features, particularly legacy function names, 
                                       // that have been replaced and are going to go away in the future.
#include <boost/filesystem/operations.hpp>
using namespace boost::filesystem;

// /* There is a better solution, if you have the Boost C++ library installed which provides us with the following function:
//  * boost::uintmax_t file_size(const path& p); 
//  */

// include some epitome headers
#include "filesStorage.hpp"
#include "makeDirectory.h"
#include "readWriteDataToFromFile.hpp"
#include "ImageEpitome.hpp"
// JPEG etc compressor to epitome reconstructed images
#include "Dis_Bpp.h"

// change row_column indices into j2k image for a less compressed file size
#include "Idx2Img.h"

// include FFT_Convolution or Correlation headers
#include "convolution_fftw.h"
#include "convolution_std.h"
// maximum in double

#define DOUBLE_MAX  DBL_MAX // no semicolons needed, since it is not a statement but a directive.
// minimum in double
#define DOUBLE_MIN  -1*DOUBLE_MAX // no semicolons needed, since it is not a statement but a directive.

using namespace std;
using namespace cv;

#endif