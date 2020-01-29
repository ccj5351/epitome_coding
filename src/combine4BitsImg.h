#ifndef _combine_4_Bits_H
#define _combine_4_Bits_H

#include <iostream>
#include <stdio.h>
#include <vector>
#include <fstream>
#include <limits.h>
#include <bitset>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp> // histogram
#include <opencv2/highgui/highgui.hpp>

using namespace std;
using namespace cv;

#if INT_MAX > 0x7FFF
typedef unsigned short uint2;  /* two-byte integer (large arrays)      */
typedef unsigned int   uint4;  /* four-byte integers (range needed)    */
#else
typedef unsigned int   uint2;
typedef unsigned long  uint4;
#endif

typedef unsigned char uchar;
/* std::size_t is the unsigned integer type of the result of the sizeof operator,
as well as the sizeof... operator and the align of operator (since C++11).
*/

// no period added at the end
#define INPUT_NUM 2
#define OUTPUT_NUM 1
#define uint unsigned int
struct combine_4_Bit_Img {
	uchar * input; //  4-valid-bit numbers
	uchar *  output; // 8-valid-bit numbers

	void initil();

	void doComine();

	void undoCombine();
	Mat_<uchar> combine_4_Bit_Img::Combine(const Mat_<uchar> & m_in);
	Mat_<uchar> combine_4_Bit_Img::UnCombine(const Mat_<uchar> & m_in);
	void clear();

};

#endif