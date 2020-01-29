#include "factorize.h"

// Code adapted from gsl/fft/factorize.c
void factorize(const int signal_length, // as an input integer to be factorized
	int * n_factors, // the number of factors
	int factors[], // the factors we get in the end
	int * implemented_factors // e.g., {13, 11, 7, 5, 3, 2, 0 } , to end with zero to detect the end of the array

	// /* In FFTW3, an optimal size is of the form 
	//  * 2^a * 3^b * 5^c * 7^d * 11^e * 13^f , with (e + f) equals 0 or 1 (see FFTW website)
	//  * We also discard the sizes multiple of 4*4*4*2 = 128, which appeared to decrease the performances.
	//  */
	){


	if (signal_length == 0)
	{
		printf("Signal length n must be positive integer!\n");
		return;
	}

	if (signal_length == 1)
	{
		factors[0] = 1;
		*n_factors = 1;
		return;
	}

	// /* deal with the implemented factors */
	int test_length = signal_length;
	int factor;
	int idx_implemented_factors = 0; // indices of the elements in the array - implemented_factors
	int idx_factors = 0; // indices of the elements in the array - factors - we get in the end
	// e.g., implemented factors = {13, 11, 7, 5, 3, 2, 0 } , to end with zero to detect the end of the array
	while (implemented_factors[idx_implemented_factors] && test_length != 1){
		factor = implemented_factors[idx_implemented_factors]; // e.g., when factor is 13, for example

		while ((test_length % factor) == 0){ // do when being divided exactly by 13, for example
			test_length = test_length / factor;
			factors[idx_factors] = factor;
			idx_factors++;
		}
		idx_implemented_factors++;
	}

	// Ok that's it
	if (test_length != 1){
		factors[idx_factors] = test_length;
		idx_factors++; 
		// The last post-increment operator is needed, for example, an array A[5], has 5 elements in total, 
		// but the index of the last element is 4, due to 0-based indices.
		// Thus, if we want to get the number of "5", we have to increment "4" by 1.
	}

	// /* check that the factorization is correct */
	{
		int product = 1;

		for (int idx = 0; idx < idx_factors; ++idx)
		{
			product *= factors[idx];
		}

		if (product != signal_length)
		{
			printf("factorization failed");
		}
	}

	*n_factors = idx_factors;
}



bool is_optimal(
	const int signal_length, // as an input integer to be factorized
	int * implemented_factors // e.g., {13, 11, 7, 5, 3, 2, 0 } , to end with zero to detect the end of the array
	){
	// to discard the sizes multiple of 4*4*4*2 = 128, which appeared to decrease the performances.
	// We check that n is not a multiple of 4*4*4*2 = 128.
	if (signal_length % 128 == 0)
		return false;

	int n_factors = 0;
	int factors[64]; // preserve a large enough array, 
	// here, this array has 64 elements, 
	// so the minimum value it can represent is 2^64 = 1.8447e+19, a very large size.

	factorize(signal_length, &n_factors, factors, implemented_factors);

	// We just have to check if the last factor belongs to FFTW_FACTORS
	// int FFTW_FACTORS[7] = {13,11,7,5,3,2,0}; // end with zero to detect the end of the array
	bool last_factor_yes = false;
	int i = 0;
	while (implemented_factors[i]) // just do, if not zero
	{
		if (factors[n_factors - 1] == implemented_factors[i])
			last_factor_yes = true;
		++i;
	}

	// if flag last_factor_yes == false, directly return false, 
	// without the following calculating about the flag num_11_13_yes
	if (!last_factor_yes)
		return false;
	else { // when the flag last_factor_yes == true

		// In FFTW, only size that can be written as 2^a * 3^b * 5^c * 7^d * 11^e * 13^f,
		// and e + f = 0 or 1, 
		// is regarded as optimal for FFT, otherwise, padding of signals is needed for FFT.
		bool num_11_13_yes = false;
		int num_11_13 = 0;
		for (int idx = 0; idx < n_factors; ++idx){
			if (factors[0] == 11 || factors[0] == 13)
				num_11_13++;	
		}

		if (num_11_13 == 0 || num_11_13 == 1)
			num_11_13_yes = true;
		
		return num_11_13_yes;
	}
}


int find_closest_factor(const int n, int * implemented_factor)
{
	if (is_optimal(n, implemented_factor))
		return n;
	else
	{
		int j = n + 1;
		while (!is_optimal(j, implemented_factor)) // do when not optimal size
			++j;
		return j;
	}
}