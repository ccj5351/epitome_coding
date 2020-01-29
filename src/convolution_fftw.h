// from https://github.com/jeremyfix/FFTConvolution

#ifndef CONVOLUTION_FFTW_H
#define CONVOLUTION_FFTW_H

#include <cassert> 
// /* C Diagnostics Library
//  * assert.h defines one macro function (i.e., assert) that can be used as a standard debugging tool :
//  * void assert (int expression);
//  * Expression to be evaluated.If this expression evaluates to 0, this causes an assertion failure that terminates the program.
//  */

#include <fftw3.h>
#include <algorithm>
#include <ctime> // show time needed for speed testing
#include <cstdlib> // (stdlib.h)
#include <string>
#include <cstring> // <cstring> (string.h). This header file defines several functions to manipulate C strings and arrays.
#include <fstream> 
#include "factorize.h"

namespace FFTW_Convolution {
	
	// /* from http://www.fftw.org/doc/Real_002ddata-DFTs.html
	//  * FFTW is best at handling sizes of the form 2a 3b 5c 7d 11e 13f,where e+f is either 0 or 1, and the other exponents are arbitrary. 
	//  * Other sizes are computed by means of a slow, general-purpose algorithm (which nevertheless retains O(n log n) performance even for prime sizes). 
	//  * (It is possible to customize FFTW for different array sizes; see Installation and Customization.) 
	//  * Transforms whose sizes are powers of 2 are especially fast, and it is generally beneficial for the last dimension of an r2c/c2r transform to be even.
	//  */
	typedef enum
	{
		// Here if no explicit "UNPADDED" occurs, we always mean PADDING as default.
		// i.e., here LINEAR_FULL == LINEAR_FULL_PADDED
		
		// *******************
		// Padding
		// *******************
		// Linear Convolution
		LINEAR_CONVOLUTION_FULL, // 0
		LINEAR_CONVOLUTION_SAME, // 1
		LINEAR_CONVOLUTION_VALID, // 2
		LINEAR_CONVOLUTION_VALIDSAME, // 3 , without circularly wrap
		
		
		// Linear Correlation
		LINEAR_CORRELATION_FULL, // 4		
		LINEAR_CORRELATION_SAME, // 5		
		LINEAR_CORRELATION_VALID, // 6		
		LINEAR_CORRELATION_VALIDSAME, // 7 , without circularly wrap
		

		// Circular Convolution
		CIRCULAR_CONVOLUTION_FULL, // 8		
		CIRCULAR_CONVOLUTION_SAME, // 9

		// *******************
		// NO Padding
		// *******************
		// Linear Convolution
		LINEAR_CONVOLUTION_FULL_UNPADDED, // 10
		LINEAR_CONVOLUTION_SAME_UNPADDED, // 11		
		LINEAR_CONVOLUTION_VALID_UNPADDED,	 // 12
		LINEAR_CONVOLUTION_VALIDSAME_UNPADDED, // 13 , without circularly wrap

		// Linear Correlation		
		LINEAR_CORRELATION_FULL_UNPADDED,	 // 14	
		LINEAR_CORRELATION_SAME_UNPADDED, // 15		
		LINEAR_CORRELATION_VALID_UNPADDED,	 // 16	
		LINEAR_CORRELATION_VALIDSAME_UNPADDED, // 17, without circularly wrap

		// Circular Convolution		
		CIRCULAR_CONVOLUTION_FULL_UNPADDED,	 // 18	
		CIRCULAR_CONVOLUTION_SAME_UNPADDED, // 19

		// 2 new kinds are added, which are only valid for the type of "VALIDSAME"
		// we bring the two new kinds, due to the circularly wrap  to
		// the convolution or correlation of "VALIDSAME"
		LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP, // 20, with circularly wrap
		LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP // 21, with circularly wrap

	} FFT_Convolution_Mode; // 20 kinds of modes in total, every item will be identified by an integer from 0 to 19, when the program is run.


	// int FFTW_FACTORS[7]; // e.g.,  = { 13, 11, 7, 5, 3, 2, 0 }; // end with zero to detect the end of the array
	// unsigned int fftw_flag; // e.g.  = FFTW_ESTIMATE; 
	// /* In FFTW3, an optimal size is of the form 
	//  * 2^a * 3^b * 5^c * 7^d * 11^e * 13^f , with (e + f) equals 0 or 1 (see FFTW website)
	//  * We also discard the sizes multiple of 4*4*4*2 = 128, which appeared to decrease the performances.
	//  */

	typedef struct FFT_Workspace
	{
		double * in_src, *in_kernel;
		double	* out_src, *out_kernel;
		int h_src, w_src;
		int h_kernel, w_kernel;
		int w_fftw, h_fftw;
		bool IsPeriodicPadding2Src;
		bool showSeconds;
	//	bool IsCircularlyWrapForValidSame;
		unsigned int fftw_flag; // e.g.  = FFTW_ESTIMATE;
		int * FFTW_FACTORS; // e.g., = { 13, 11, 7, 5, 3, 2, 0 }; // end with zero to detect the end of the array
		FFT_Convolution_Mode mode;
		std::string s_mode;
		double * dst_ifft; // normalized IFFT, i.e., normalized raw convolution, without indices selection
		double * dst_convolution; // convolution corresponding to CONVOLUTION MODE, i.e., the selected of normalized "dst_ifft"
		int h_dst, w_dst; // the size of convolution output ; This is automatically set by init_workspace
		fftw_plan p_forw_src;
		fftw_plan p_forw_kernel;
		fftw_plan p_back;
	} FFT_Workspace; // do not worry the Workspace in FFTW has the same name as that in standard convolution
	// for the using of namespace FFTW_Convolution, and namespace STD_Convolution.

	void init_workspace(FFT_Workspace & ws, FFT_Convolution_Mode mode, int h_src, int w_src, int h_kernel, int w_kernel, bool verbose);
	void update_workspace_few(FFT_Workspace & ws);

	void clear_workspace(FFT_Workspace & ws);

	// Only do FFT to src, instead of convolution or correlation
	void fftw_src_fft(FFT_Workspace &ws, double * src);
	// only do FFT to Kernel, instead of convolution or correlation
	void fftw_kernel_fft(FFT_Workspace &ws, double * kernel);

	// to calculate convolution or correlation from two FFT signal, not the original time-domain or space domain signals
	void fftw_circular_convolution_from_FFT_Signal(FFT_Workspace &ws);

	// Before calling this function, we must have gotten the FFT signals of src and kernel,
	// i.e., the FFT signal of src has been saved into ws.out_src, 
	// and the FFT signal of kernel has been saved into ws.out_kernel;
	// Otherwise, we will encounter errors !!!
	void convolve_from_FFT_Signal(FFT_Workspace &ws);

	// Compute the circular convolution of src and kernel modulo ws.h_fftw, ws.w_fftw
	// using the Fast Fourier Transform (FFT)
	// The result is in ws.dst
	void fftw_circular_convolution(FFT_Workspace &ws, double * src, double * kernel);

	void convolve(FFT_Workspace &ws, double * src, double * kernel);

	// /*  The relationship between correlation and convolution
	//  *  let the MATLAB function of xcorr(a, b) be linear correlation of a and b,
	//  *  conv(a,b) be linear convolution of a and b,
	//  *  cconv(a,b) be circular convolution of a and , respectively.
	//  *   1) calculates circular convolution:
	//  *        * do circular convolution, via some length L_fft fft and ifft, i.e., ifft{[fft(src,L_fft).*fft(kernel,L_fft)], L_fft}
	//  *   2) calculates linear convolution via circular convolution:
	//  *        * when the FFT/IFFT length L_fft is equal to or larger than L = length(src) + length(kernel) -1, 
	//  *        * the circular convolution result then is the same as the linear convolution result that we want to get before.
	//  *   3) calculates linear convolution via linear correlation 
	//  *      NOTE: AND you should know that the linear convolution is always to be got via CIRCULAR convolution at last):
	//  *        * take 1-D or 2-D real array for example,
	//  *        * 1-D array: xcorr(a1[n], b1[n]) = conv(a1[n], b1[-n]) = conv(a1[n], flip(b1[n]));
	//  *        * 2-D array: xcorr(a2[n,m], b2[n,m]) = conv(a2[n,m], b2[-n,-m]) = conv(a2[n,m], flip(flip(b2[n,m],dim-column), dim-row));
	//  *        * Thus, if we want to get correlation, first we completely flip the kernel - b1[n] or b2[n,m] for example here, 
	//  *        * then call the linear convolution function to the input signal and the flipped kernel.
	//  *        * Definitely, the so-called linear convolution will have to be calculated via the above CIRCULAR convolution at last.


	void displayConvolve(FFT_Workspace & ws);

	void saveConvolve(FFT_Workspace & ws, std::string & filename);

	void help();

}
#endif
