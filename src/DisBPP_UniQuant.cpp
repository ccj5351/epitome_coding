#include "Dis_Bpp.h"

int Distortion_Bpp::UniformQuantize(const Mat_<double> & src, Mat_<unsigned int> & dst, UniformQuant & uniQuant){

	if (!src.data){
		cerr << "No Valid Input Image Being Read Now!\n";
		return 1;
	}
	else{
		for (int i = 0, i_size = src.rows; i != i_size; ++i){
			for (int j = 0, j_size = src.cols; j != j_size; ++j){
				dst.at<unsigned int>(i, j) = uniQuant.Symbol(src.at<double>(i, j));
			}
		}

		return 0;
	}
}

int Distortion_Bpp::deUniformQuantize(const Mat_<unsigned int> & src, Mat_<double> & dst, UniformQuant & uniQuant){
	if (!src.data){
		cerr << "No Valid Input Image Being Read Now!\n";
		return 1;
	}
	else{
		for (int i = 0, i_size = src.rows; i != i_size; ++i){
			for (int j = 0, j_size = src.cols; j != j_size; ++j){
				dst.at<double>(i, j) = uniQuant.Value(src.at<unsigned int>(i, j));
			}
		}
		return 0;
	}
}

int Distortion_Bpp::UniformQuantize(const Mat_<double> & src, Mat_<BYTE> & dst, UniformQuant & uniQuant){
	if (!src.data){
		cerr << "No Valid Input Image Being Read Now!\n";
		return 1;
	}
	else{
		for (int i = 0, i_size = src.rows; i != i_size; ++i){
			for (int j = 0, j_size = src.cols; j != j_size; ++j){
				dst.at<uchar>(i, j) = uniQuant.Symbol(src.at<double>(i, j));
			}
		}

		return 0;
	}
}

int Distortion_Bpp::deUniformQuantize(const Mat_<BYTE> & src, Mat_<double> & dst, UniformQuant & uniQuant){
	if (!src.data){
		cerr << "No Valid Input Image Being Read Now!\n";
		return 1;
	}
	else{
		for (int i = 0, i_size = src.rows; i != i_size; ++i){
			for (int j = 0, j_size = src.cols; j != j_size; ++j){
				dst.at<double>(i, j) = uniQuant.Value(src.at<uchar>(i, j));
			}
		}
		return 0;
	}
}

