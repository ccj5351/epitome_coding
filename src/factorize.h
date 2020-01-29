// from https://github.com/jeremyfix/FFTConvolution

#ifndef FACTORIZE_H
#define FACTORIZE_H

#include <iostream>

// Code adapted from gsl/fft/factorize.c
void factorize(const int signal_length, // as an input integer to be factorized
	int * n_factors, // the number of factors
	int factors[], // the factors we get in the end
	int * implemented_factors // e.g., {13, 11, 7, 5, 3, 2, 0 } , to end with zero to detect the end of the array

	// /* In FFTW3, an optimal size is of the form 
	//  * 2^a * 3^b * 5^c * 7^d * 11^e * 13^f , with (e + f) equals 0 or 1 (see FFTW website)
	//  * We also discard the sizes multiple of 4*4*4*2 = 128, which appeared to decrease the performances.
	//  */
	);



bool is_optimal(
	const int n, // as an input integer to be factorized
	int * implemented_factors // e.g., {13, 11, 7, 5, 3, 2, 0 } , to end with zero to detect the end of the array
	);


int find_closest_factor(const int n, int * implemented_factor);

#endif
