/*---------------------------------------------------------------------------*/
// Image Compression Toolbox v1.2
// written by
// Satish Kumar S
// satishkumr@lycos.com
//
// Copyright 1999 Satish Kumar S
//
// Permission is granted to use this software for research purposes as
// long as this notice stays attached to this software.
/*---------------------------------------------------------------------------*/

#ifndef _QUANTIZER_H
#define _QUANTIZER_H

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <stdlib.h>
#include "bitstream.h"
#include "opencv2/core/core.hpp"
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
using namespace std;
using namespace cv;

#define uint unsigned int

class UniformQuant{

public:
    double minvalue, stepsize, invstepsize;
    int nLevels;

public:
	UniformQuant();

	UniformQuant(double maxv, double minv, int nlev);

	UniformQuant(double min, double stsize);
	

	// encode a real value into its quantized value, represented by one code or symbol of the codebook or dictionary.
    inline unsigned int Symbol(double value){
        return ((unsigned int)((value - minvalue) * invstepsize));
    }

	// decode a symbol into its original approximate real value, represented by the median of the interval that the original real value falls in.
    inline double Value(unsigned int code){
        return (minvalue + ((double)code + 0.5)*stepsize); // "0.5" makes we get the median of the interval;
    }

	// clear existing parameters, and set them being zeros
	inline void clearParams(){
		minvalue = stepsize = invstepsize = 0.0;
		nLevels = 0;
	}

	// display parameters
	inline void displayParams(){
		std::cout << "minvalue = " << minvalue
			<< ", stepsize = " << stepsize
			<< ", nLevels = " << nLevels << std::endl;
	}

	// write to file the important value, which are necessary for the following de-quantization
	void WriteHeader(BitStream & bs);

	// read the header for loading the important values for the following de-quantization
	void ReadHeader(BitStream &bs);
	// do uniform quantization
	void UniformQuantize(const Mat_<double> & src, Mat_<uchar> & dst);
	void UniformQuantize(const std::vector<int> & src, std::vector<int> & dst);
	// do uniform de-quantization
	void deUniformQuantize(const Mat_<uchar> & src, Mat_<double> & dst);
	void deUniformQuantize(const std::vector<int> & src, std::vector<double> & dst);
	void deUniformQuantize(const std::vector<int> & src, std::vector<int> & dst);
};

#endif