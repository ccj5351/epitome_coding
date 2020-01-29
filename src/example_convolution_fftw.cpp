#define _FFTW_CONVOLUTION_EXAMPLE_
#ifndef  _FFTW_CONVOLUTION_EXAMPLE_

#include "convolution_fftw.h"
#include <ctime>
#include <fstream>
using namespace std;
using namespace FFTW_Convolution;


//***************************************
//**********  main functioin  ***********
//***************************************
int main(int argc, char *argv[]){
	// clock_t clock (void), Returns the processor time consumed by the program.
	// The value returned is expressed in clock ticks, which are units of time of a constant but system - specific length(with a relation of CLOCKS_PER_SEC clock ticks per second).
	// The epoch used as reference by clock varies between systems, but it is related to the program execution(generally its launch).
	// To calculate the actual processing time of a program, the value returned by clock shall be compared to a value returned by a previous call to the same function.
	
	
	std::cout << " Calculating...\n";
	clock_t t = clock();
	FFT_Workspace ws;
	FFT_Convolution_Mode mode; // e.g. = LINEAR_CONVOLUTION_SAME;
	
	int IterationNum = 1;
	int h_src = 1;
	int w_src = 5;
	int w_kernel = 3;
	int h_kernel = 1;
	bool verbose = true;

	help();

	for (int i = 0; i != 20; ++i){
		mode = (FFT_Convolution_Mode)i;
	init_workspace(ws, mode, h_src,w_src, h_kernel, w_kernel, verbose);
	string filename = "E:/imageDatabase/fft_result/" + ws.s_mode + ".txt";
	double  src[] = { 1, 2, 3, 4, 5 };
	double  kernel[] = { 1, 2, 3 };
	
	// method 1, to calculate convolution or correlation
	for (int i = 0; i < IterationNum; ++i){
	convolve(ws, src, kernel);
	displayConvolve(ws);
	//saveConvolve(ws, filename);
	}

	std::cout << "\n\n";

	double  src1[] = { 1, 2, 3, 4, 5 };
	double  kernel1[] = { 1, 2, 3 };
	string filename1 = "E:/imageDatabase/fft_result/" + ws.s_mode + "1.txt";
	for (int i = 0; i < IterationNum; ++i){

		// method 2, to calculate convolution or correlation
		fftw_src_fft(ws, src1);
		fftw_kernel_fft(ws, kernel1);
		convolve_from_FFT_Signal(ws);
		displayConvolve(ws);
		//saveConvolve(ws, filename1);
	}

	std::cout << "\n\n\n\n";

	}

	t = clock() - t;
	printf("It took me %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);

	clear_workspace(ws);

	return 0;
}
#endif