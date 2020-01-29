// from https://github.com/jeremyfix/FFTConvolution

#ifndef CONVOLUTION_STD_H
#define CONVOLUTION_STD_H
#include <algorithm> // e.g., using std::max, std::min, etc
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>


// /*
//  * #define max(a,b)            (((a) > (b)) ? (a) : (b))
//  * #define min(a,b)            (((a) < (b)) ? (a) : (b))
//  */

namespace STD_Convolution
{

	enum Std_Convolution_Mode
	{   // Linear Convolution
		LINEAR_CONVOLUTION_FULL, // 0
		LINEAR_CONVOLUTION_SAME, // 1
		LINEAR_CONVOLUTION_VALID, // 2
		LINEAR_CONVOLUTION_VALIDSAME, // 3 

		// Linear Correlation
		LINEAR_CORRELATION_FULL, // 4		
		LINEAR_CORRELATION_SAME, // 5		
		LINEAR_CORRELATION_VALID, // 6		
		LINEAR_CORRELATION_VALIDSAME, // 7

		// Circular Convolution
		CIRCULAR_CONVOLUTION_FULL, // 8		
		CIRCULAR_CONVOLUTION_SAME, // 9
	} ;

	 struct Std_Workspace
	{
		int h_src, w_src; // height and width of the input image
		int h_kernel, w_kernel; // height and width of the kernel

		int offset_h_dst, offset_w_dst; // shifting variable compared with the FULL case
		// for full case, they should be zeros;
		// for same case, they should be int(ws.h_kernel / 2.0) and int(ws.w_kernel / 2.0), respectively;
		// for valid case, they should be (ws.h_kernel - 1) and(ws.w_kernel - 1), respectively;
		// for VALIDSAME case, they should be (ws.h_kernel - 1) and(ws.w_kernel - 1), respectively, the same as to the valid case,
		// but the h_dst and w_dst are different from that of the case of VALID.
		Std_Convolution_Mode mode;
		bool showSeconds;
		std::string s_mode;
		double * dst;
		int h_dst, w_dst; // height and width of the output image
	} ;

	 void init_workspace(Std_Workspace & ws, Std_Convolution_Mode mode, int h_src, int w_src, int h_kernel, int w_kernel, bool verbose);
	

	 void clear_workspace(Std_Workspace & ws); // to delete output array


	 void update_workspace(Std_Workspace & ws, Std_Convolution_Mode mode, int h_src, int w_src, int h_kernel, int w_kernel);
	 
	 void general_convolve(Std_Workspace & ws, double * src, double * kernel);

	 void convolve(Std_Workspace &ws, double * src, double * kernel);
	 void saveConvolve(Std_Workspace & ws, std::string & filename);

	

} // end of namespace STD_Convolution


#endif
