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

#include "quantizer.h"

using namespace std;
using namespace cv;

UniformQuant::UniformQuant(){
		minvalue = stepsize = invstepsize = 0.0;
		nLevels = 0;
	}


UniformQuant::UniformQuant(double maxv, double minv, int nlev){

		minvalue = minv;
		nLevels = nlev;

		if (nlev == 1){            // only a single level??        
			stepsize = (maxv - minv);
			invstepsize = 0.0;
		}

		else{
			stepsize = (maxv - minv) / (nlev - 1);
			invstepsize = 1.0 / stepsize;
		}
	}

UniformQuant::UniformQuant(double min, double stsize){
		minvalue = min;
		stepsize = stsize;
		invstepsize = 1.0 / stepsize;
	}


	// write to file the important value, which are necessary for the following de-quantization
void UniformQuant::WriteHeader(BitStream & bs){
		// The size of an int is really implementation dependent. 
		// Back in the day, when processors were 16 bit, an int was 2 bytes. Nowadays, it's most often 4 bytes (32 bits).
		// Still, using sizeof(int) is the best way to get the size of an integer.
		size_t int_size = sizeof(int) * 8; // e.g., int_size = 32 bits = 4 bytes;
#ifdef _DEBUG
		cout << "size of int = " << int_size << "bits.\n";
#endif
		bs.WriteBits(nLevels, int_size);
		if (nLevels <= 1){
			float meanvalue = (float)(minvalue + 0.5*stepsize);
			bs.WriteBits(*((unsigned int *)&meanvalue), int_size);
		}

		else{
			float minv = (float)minvalue;
			float stsize = (float)stepsize;
			bs.WriteBits(*((unsigned int *)&minv), int_size);
			bs.WriteBits(*((unsigned int *)&stsize), int_size);
		}
	}

	// read the header for loading the important values for the following de-quantization
void UniformQuant::ReadHeader(BitStream &bs){
		size_t int_size = sizeof(int) * 8; // e.g., int_size = 32;
#ifdef _DEBUG
		cout << "size of int = " << int_size << "bits.\n";
#endif
		nLevels = bs.ReadBits(int_size);
		if (nLevels == 1){
			unsigned int tmeanv = bs.ReadBits(int_size);
			minvalue = *(float *)&tmeanv;
			stepsize = 0.0;
			invstepsize = 0.0;
		}
		else{
			unsigned int tminv = bs.ReadBits(int_size);
			unsigned int tsteps = bs.ReadBits(int_size);
			minvalue = *(float *)&tminv;
			stepsize = *(float *)&tsteps;
			invstepsize = 1.0 / stepsize;
		}
}


void UniformQuant::UniformQuantize(const Mat_<double> & src, Mat_<uchar> & dst){
		if (!src.data){
			cerr << "No Valid Input Image Being Read Now!\n";
		}
		else{
			for (int i = 0, i_size = src.rows; i != i_size; ++i){
				for (int j = 0, j_size = src.cols; j != j_size; ++j){
					dst.at<uchar>(i, j) = Symbol(src.at<double>(i, j));
				}
			}

		}
	}

// overload + 1
void UniformQuant::UniformQuantize(const std::vector<int> & src, std::vector<int> & dst){
	if (src.empty()){
		cerr << "No Valid Input Image Being Read Now!\n";
	}
	else{
		for (int i = 0, i_size = src.size(); i != i_size; ++i){
				dst[i] = Symbol(double(src[i]));
		}
	}
}

void UniformQuant::deUniformQuantize(const Mat_<uchar> & src, Mat_<double> & dst){
		if (!src.data){
			cerr << "No Valid Input Image Being Read Now!\n";
		}
		else{
			for (int i = 0, i_size = src.rows; i != i_size; ++i){
				for (int j = 0, j_size = src.cols; j != j_size; ++j){
					dst.at<double>(i, j) = Value(src.at<uchar>(i, j));
				}
			}
		}
}

// overload + 1
void UniformQuant::deUniformQuantize(const std::vector<int> & src, std::vector<double> & dst){
	if (src.empty()){
		cerr << "No Valid Input Image Being Read Now!\n";
	}
	else{
		for (int i = 0, i_size = src.size(); i != i_size; ++i){
				dst[i] = Value(src[i]);	
		}
	}
}
// overload + 2
void UniformQuant::deUniformQuantize(const std::vector<int> & src, std::vector<int> & dst){
	if (src.empty()){
		cerr << "No Valid Input Image Being Read Now!\n";
	}
	else{
		for (int i = 0, i_size = src.size(); i != i_size; ++i){
			dst[i] = (int)Value(src[i]);
		}
	}
}
