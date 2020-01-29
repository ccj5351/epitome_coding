#include "allHeaders.hpp"
#include "Dis_Bpp.h"
#include "quantizer.h"

void Distortion_Bpp::doBinUniformQuantCodec(
	// const int & jpeg_quality, //  quantization level belonging to [1, 100], and higher values mean higher resolution;
	const Mat_<double> & src, // original image
	const Mat_<double> & epi_recon, // epitome reconstruction
	const bool isConstQuantizationLevel, // set a quantization level or not
	const int QuantizationLevel, // constant quantization level value, 0 <= QuantizationLevel <= 100;
	const string & binary_error_save_path, // indicating binary file for saving the quantization result
	const string & final_recon_img_name, // final reconstruction result from epitome reconstruction and uniform quantization of residual images.
	const bool & IsSaveImgFlag,
	double & psnr,
	double & mse,
	std::fstream * fs
	 
	){

	char c_jpeg_quality[10];
	// set three-character length
	sprintf(c_jpeg_quality, "%03d", jpeg_quality);
	string s_jpeg_quality(c_jpeg_quality);
	Mat_<double> error = src - epi_recon;

	// keep the max and min values of error, before normalization.
	double max_intensity, min_intensity;
	cv::minMaxLoc(error, &min_intensity, &max_intensity);
#ifdef _DEBUG
	std::cout << "min = " << min_intensity << ", max = " << max_intensity << "\n";
#endif

	//********************************************
	// some variables for uniform quantization
	//********************************************

	// here we use uchar, since we will guarantee that the quantized value level is equal to or less than 256.
	Mat_<uchar> quantized_unnormalized_error = Mat::zeros(src.rows, src.cols, CV_8UC1);
	int nLevels;
	// Function : FILE * fopen ( const char * filename, const char * mode );
	FILE * p_file = NULL;
	BitStream bs(p_file);
	if (isConstQuantizationLevel){
		nLevels = (int)(QuantizationLevel * (max_intensity - min_intensity) / 100);
	}
	else{
		if (jpeg_quality == 0){
			nLevels = 1; // 1<= nLevels <= 256;
		}
		else{
			nLevels = (int)(jpeg_quality * (max_intensity - min_intensity) / 100);
		}
	}

	// at last, make sure 1 <= nLevels <= 256
	nLevels = nLevels > 256 ? 256 : nLevels;
	nLevels = nLevels < 1 ? 1 : nLevels;

	int bit_for_quantized_error;
	if (nLevels == 1)
		bit_for_quantized_error = 1;
	else 
		bit_for_quantized_error = (int)ceil(log2(nLevels));

	UniformQuant uniQuant_BINARY_FILE(max_intensity, min_intensity, nLevels);
	char * char_error_binary_file_path = new char[binary_error_save_path.length() + 1];
	std::strcpy(char_error_binary_file_path, binary_error_save_path.c_str());

	//**********************************************************
	// encoding or writing quantization result to binary files;
	//**********************************************************
	p_file = fopen(char_error_binary_file_path, "wb");

	if (!p_file){ // if NULL
		std::cout << "Cannot open the file of " << char_error_binary_file_path << endl;
	}
	else{
		bs.SetFilePointer(p_file);
		//************************************
		// for writing header binary file
		//************************************
		uniQuant_BINARY_FILE.WriteHeader(bs);

		//************************************
		// for body binary file
		//************************************
		uniQuant_BINARY_FILE.UniformQuantize(error, quantized_unnormalized_error);
		for (int i = 0, i_size = quantized_unnormalized_error.rows; i != i_size; ++i){
			for (int j = 0, j_size = quantized_unnormalized_error.cols; j != j_size; ++j){
				bs.WriteBits(quantized_unnormalized_error.at<uchar>(i, j), bit_for_quantized_error); // bit_for_quantized_error
			}
		}

#ifdef _DEBUG
		std::cout << "Before bs.Flush(), numBitsWritten = " << bs.numBitsWritten << endl;
#endif
		bs.Flush();
#ifdef _DEBUG
		std::cout << "After  bs.Flush(), numBitsWritten = " << bs.numBitsWritten << endl;
#endif

		/*
		function : fclose
		int fclose ( FILE * stream );
		Closes the file associated with the stream and disassociates it.
		All internal buffers associated with the stream are disassociated from it and flushed,
		i.e.,the content of any unwritten output buffer is written and the content of any unread input buffer is discarded.
		*/
		fclose(p_file);
	} /*end of encoding*/


	//***************************************************************************
	// decoding or reading de-quantization result from physical binary files;
	//***************************************************************************
	p_file = fopen(char_error_binary_file_path, "rb");
	if (!p_file){ // if NULL
		cout << "Cannot open the file of " << char_error_binary_file_path << endl;
	}
	else{
		cout << "Successfully open the file of " << char_error_binary_file_path << endl;
		fseek(p_file, 0, SEEK_SET);
		bs.clear();
		bs.SetFilePointer(p_file);
#ifdef _DEBUG
		cout << "Before Clearing parameters:" << endl;
		uniQuant_BINARY_FILE.displayParams();
#endif
		uniQuant_BINARY_FILE.clearParams();
#ifdef _DEBUG
		cout << "Before ReadHeader:" << endl;
		uniQuant_BINARY_FILE.displayParams();
#endif
		//************************************
		// for reading header binary file
		//************************************
		uniQuant_BINARY_FILE.ReadHeader(bs);

#ifdef _DEBUG
		cout << "After ReadHeader:" << endl;
		uniQuant_BINARY_FILE.displayParams();
#endif
		//************************************
		// for body binary file
		//************************************
		bs.clear();
		bs.SetFilePointer(p_file);
		Mat_<uchar> decode_error = Mat_<uchar>::zeros(src.rows, src.cols);
#ifdef _DEBUG
		cout << "Before Reading, numBitsRead = " << bs.numBitsRead << endl;
#endif
		for (int i = 0, i_size = decode_error.rows; i != i_size; ++i){
			for (int j = 0, j_size = decode_error.cols; j != j_size; ++j){
				decode_error.at<uchar>(i, j) = bs.ReadBits(bit_for_quantized_error);
			}
		}

		fclose(p_file);
		error = Mat_<double>::zeros(src.rows, src.cols);
		uniQuant_BINARY_FILE.deUniformQuantize(decode_error, error);
#ifdef _DEBUG
		cout << "After Reading, numBitsRead = " << bs.numBitsRead << endl;
#endif

		Mat_<double> recon_src = Mat_<double>(error + epi_recon);
		double max_intensity = .0, min_intensity = .0;
		cv::minMaxLoc(recon_src, &min_intensity, &max_intensity);
		if ((min_intensity < .0) | (max_intensity > 255.0)){/*If Intensities overflow*/
			std::cout << "Intensities overflow!\n";
			// intensity should belong to the range [0, 255];
			for (int i = 0, rowNum = recon_src.rows; i != rowNum; ++i)
				for (int j = 0, colNum = recon_src.cols; j != colNum; ++j){
				recon_src.at<double>(i, j) = std::max<double>(recon_src.at<double>(i, j), 0.0);// >= 0;
				recon_src.at<double>(i, j) = std::min<double>(recon_src.at<double>(i, j), 255.0);// <= 255
				}
		}/*End of If Intensities overflow*/

		if (IsSaveImgFlag){
		imwrite(final_recon_img_name, recon_src);
		}

		getEpiPSNR_MSE(src, recon_src, psnr, mse);
		if (fs->is_open()){
			*fs << "Uniform Quantization level = " << s_jpeg_quality << ", nlevels = " << nLevels
				<< ", bits = " << bit_for_quantized_error << ", psnr = " << psnr << ", mse = " << mse << endl;
			std::cout << "Operation successfully performed\n";
		}
		else{
			std::cout << "Error opening std::fstream file";
		}
	} /*end of decoding*/

	delete char_error_binary_file_path;
}