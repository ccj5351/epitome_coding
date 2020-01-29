#define _example_convolve_time_test_
#ifndef  _example_convolve_time_test_

#include "convolution_fftw.h"
#include "convolution_std.h"
#include <ctime>
#include <fstream>
#include <opencv2/core/core.hpp>

#define VERBOSE false
#define SAVE_RESULTS true

using namespace std;
namespace c_fft = FFTW_Convolution;
namespace c_std = STD_Convolution;


//***************************************
//**********  main functioin  ***********
//***************************************
int main(int argc, char *argv[]){
	// clock_t clock (void), Returns the processor time consumed by the program.
	// The value returned is expressed in clock ticks, which are units of time of a constant but system - specific length(with a relation of CLOCKS_PER_SEC clock ticks per second).
	// The epoch used as reference by clock varies between systems, but it is related to the program execution(generally its launch).
	// To calculate the actual processing time of a program, the value returned by clock shall be compared to a value returned by a previous call to the same function.


	std::cout << " Calculating...\n";

	//*************************************
	// signals of src and kernel
	//*************************************
	// NOTE: 
	// before we call convole() funtion, please redefine the signal src and kernel
	// for example, if the some function wants to calculate the correlation, then the kernel might have been reversed or flipped,
	// thus, for a new calculation, we have to redefine the signal, if case of some unexpected change occurred to src or kernel.
	int IterationTimes = 10000; // for iterations
	int h_src, w_src, h_kernel, w_kernel;
	w_src = atoi(argv[1]);
	h_src = w_src;
	w_kernel = atoi(argv[2]);
	h_kernel = w_kernel;


	// Build random images to convolve
	double * src = new double[h_src*w_src];
	double * kernel = new double[h_kernel*w_kernel];
	for (int i = 0; i < h_src; ++i)
		for (int j = 0; j < w_src; ++j)
			src[i*w_src + j] = rand() / double(RAND_MAX);
	// build a 2-d symmetric kernel
	// so flip the kernel, the result is the same.
	for (int i = 0; i < h_kernel; ++i)
		for (int j = 0; j < w_kernel; ++j)
			kernel[i*w_kernel + j] = 1.0;
	
	if (VERBOSE) printf("Image of size %i x %i , kernel of size %i x %i \n", h_src, w_src, h_kernel, w_kernel);

	//*************************************
	// for fftw_convolution
	//*************************************
	clock_t t = clock();
	c_fft::FFT_Workspace fft_ws;
	c_fft::FFT_Convolution_Mode fft_mode = (c_fft::FFT_Convolution_Mode)atoi(argv[3]);
	bool verbose = false;
	c_fft::init_workspace(fft_ws, fft_mode, h_src, w_src, h_kernel, w_kernel, verbose);
	string filename = "E:/imageDatabase/fft_result/" + fft_ws.s_mode + ".txt";
	if (VERBOSE) printf("Now calculating convolution via FFTW\n");
	for (int i = 0; i < IterationTimes; ++i){
		c_fft::convolve(fft_ws, src, kernel);
		// c_fft::displayConvolve(fft_ws);
		
	}
	t = clock() - t;
	float t_fftw = ((float)t) / CLOCKS_PER_SEC;
	if (VERBOSE) printf("FFTW_Convolve took me %f seconds.\n", t_fftw);
	

	//*************************************
	// for standard_convolution
	//*************************************

	t = clock() - t;
	c_std::Std_Workspace std_ws;
	c_std::Std_Convolution_Mode std_mode;
	// And compute the linear convolution
	std_mode = (c_std::Std_Convolution_Mode)atoi(argv[3]);
	c_std::init_workspace(std_ws, std_mode, h_src, w_src, h_kernel, w_kernel, verbose);
	for (int i = 0; i < IterationTimes; ++i){
		c_std::convolve(std_ws, src, kernel);
	}
	t = clock() - t;
	float t_std = ((float)t) / CLOCKS_PER_SEC;
	if (VERBOSE) printf("Standard_Convolve took %f seconds.\n", t_std);

	
	if (SAVE_RESULTS){
		printf(" The results will be saved in %s \n", filename);
		c_fft::saveConvolve(fft_ws, filename);
		c_std::saveConvolve(std_ws, filename);
	}
	
	string filename1 = "E:/imageDatabase/fft_result/timeResults.txt";
	std::ofstream fout(filename1, std::ios::app);
	if (!fout){
		std::cout << "File Not Opened" << std::endl;
	}
	fout << "Standard_Convolve took " << t_std << " seconds.\n"
		<< "FFTW_Convolve took " << t_fftw << " seconds.\n";
	
	delete[] src;
	delete[] kernel;
	c_fft::clear_workspace(fft_ws);
	c_std::clear_workspace(std_ws);
	return 0;
}
#endif