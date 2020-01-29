#include "combine4BitsImg.h"

// combine_4_Bits

void combine_4_Bit_Img::initil(){
	input = new uchar[INPUT_NUM];
	output = new uchar[OUTPUT_NUM];
}

void combine_4_Bit_Img::doComine(){
	output[0] = ((input[1] & 0xF0) >> 4) | (input[0] & 0xF0); // the high 4 bits of input[0] and the high 4 bits of input[1]
}


void combine_4_Bit_Img::undoCombine(){
	input[0] = (output[0] & 0xF0);
	input[1] = (output[0] & 0x0F) << 4;
}

Mat_<uchar> combine_4_Bit_Img::Combine(const Mat_<uchar> & m_in){
	uint height = m_in.rows,
		width = m_in.cols;
	if (width % 2 != 0)
		width += 1;
	Mat_<uchar> m_out = Mat_<uchar>::zeros(height, width / 2);
	for (int i = 0; i < height; ++i){
		for (int j = 0, k = 0; j < width; j += 2, ++k){
			input[0] = m_in.at<uchar>(i, j);
			input[1] = m_in.at<uchar>(i, j + 1);
			doComine();
			m_out.at<uchar>(i, k) = output[0];
		}
	}
	return m_out;
}

Mat_<uchar> combine_4_Bit_Img::UnCombine(const Mat_<uchar> & m_in){
	int height = m_in.rows,
		width = 2* m_in.cols;
	Mat_<uchar> m_out = Mat_<uchar>::zeros(height, width);
	for (int i = 0; i < height; ++i){
		for (int j = 0, size_j = m_in.cols; j < size_j; ++j){
			output[0] = m_in.at<uchar>(i, j);
			undoCombine();
			m_out.at<uchar>(i, 2*j) = input[0];
			m_out.at<uchar>(i, 2*j + 1) = input[1];
		}
	}
	return m_out;
}



void combine_4_Bit_Img::clear(){
	delete[] input;
	delete[] output;
}