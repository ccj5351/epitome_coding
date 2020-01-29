#include "allHeaders.hpp"
#include "Dis_Bpp.h"
#include "quantizer.h"

// Usage :
// Through uniform quantization, this function will encode the residual of src and epi_recon,
// into the result named as error_img_name(i.e., residual image) 
// and binary_error_save_path(i.e., quantization parameters necessary for the following decoding).
void Distortion_Bpp::doBMPUniformQuantEncoder(
	// const int & jpeg_quality, //  quantization level belonging to [1, 100], and higher values mean higher resolution;
	const Mat_<double> & src, // original image
	const Mat_<double> & epi_recon, // epitome reconstruction
	const bool isConstQuantizationLevel, // set a quantization level or not
	const int ConstQuantizationLevel, // constant quantization level value, 0 <= QuantizationLevel <= 100;
	const string & binary_error_save_path, // indicating binary file for saving the quantization result
	const string & error_img_name //  Uniform quantization result of residual images in the BMP format.
	){

	char * c_jpeg_quality = new char[10];
	// char c_jpeg_quality[10];
	// set three-character length
	sprintf(c_jpeg_quality, "%03d", jpeg_quality);
	string s_jpeg_quality(c_jpeg_quality);
	Mat_<double> error(src - epi_recon);

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
		nLevels = (int)(ConstQuantizationLevel * (max_intensity - min_intensity) / 100.0);
	}
	else{
		if (jpeg_quality == 0){
			nLevels = 1; // 1<= nLevels <= 256;
		}
		else{
			nLevels = (int)(jpeg_quality * (max_intensity - min_intensity) / 100.0);
		}
	}

	// at last, make sure 1 <= nLevels <= 256
	nLevels = nLevels > 256 ? 256 : nLevels;
	nLevels = nLevels < 1 ? 1 : nLevels;

	UniformQuant uniQuant_FILE(max_intensity, min_intensity, nLevels);
	// function: std::string::c_str
	// syntax: const char* c_str() const noexcept;
	// Returns a pointer to an array that contains a null-terminated sequence of characters (i.e., a C-string) 
	// representing the current value of the string object.
	// The pointer returned may be invalidated by further calls to other member functions that modify the object.

	char * char_error_binary_file_path = new char [binary_error_save_path.length() + 1];
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
		uniQuant_FILE.WriteHeader(bs);

		/*
		function : fclose
		int fclose ( FILE * stream );
		Closes the file associated with the stream and disassociates it.
		All internal buffers associated with the stream are disassociated from it and flushed,
		i.e.,the content of any unwritten output buffer is written and the content of any unread input buffer is discarded.
		*/
		fclose(p_file);

		//************************************
		// No body binary file, instead, to generate an ".bmp" image.
		//************************************
		uniQuant_FILE.UniformQuantize(error, quantized_unnormalized_error);
		// make sure the image in the format of BMP
		string temp = error_img_name;
		setFileExtension(temp, ".bmp");
		imwrite(temp, quantized_unnormalized_error);

	} /*end of encoding*/

	// memory release
	delete [] c_jpeg_quality;
	delete [] char_error_binary_file_path;
	
}


void Distortion_Bpp::doBMPUniformQuantDecoder(
	// const int & jpeg_quality, //  quantization level belonging to [1, 100], and higher values mean higher resolution;
	const Mat_<double> & src, // original image
	const Mat_<double> & epi_recon, // epitome reconstruction
	const string & binary_error_save_path, // indicating binary file for saving the quantization result
	const string & error_img_name, //  Uniform quantization result of residual images in the BMP format.
	const string & recon_img_name, // final reconstruction result from epitome reconstruction and uniform quantization of residual images.
	const bool & IsSaveImgFlag,
	double & psnr,
	double & mse,
	std::fstream * fs
	){

	//********************************************
	// some variables for uniform quantization
	//********************************************
	char * c_jpeg_quality = new char [10];
	// set three-character length
	sprintf(c_jpeg_quality, "%03d", jpeg_quality);
	string s_jpeg_quality(c_jpeg_quality);

	// here we use uchar, since we will guarantee that the quantized value level is equal to or less than 256.
	Mat_<uchar> quantized_error = imread(error_img_name, img_read_flag);
	// Function : FILE * fopen ( const char * filename, const char * mode );
	FILE * p_file = NULL;
	BitStream bs(p_file);
	UniformQuant uniQuant_BINARY_FILE;

	char * char_error_binary_file_path = new char[binary_error_save_path.length() + 1];
	std::strcpy(char_error_binary_file_path, binary_error_save_path.c_str());

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
		
		int bit_for_quantized_error;
		if (uniQuant_BINARY_FILE.nLevels == 1)
			bit_for_quantized_error = 1;
		else
			bit_for_quantized_error = (int)ceil(log2(uniQuant_BINARY_FILE.nLevels));


		//************************************
		// for body
		//************************************
		Mat_<double> decoded_error(src.rows, src.cols);
		uniQuant_BINARY_FILE.deUniformQuantize(quantized_error, decoded_error);
		
		Mat_<double> recon_src(src.rows, src.cols); 
		// = Mat_<double>::zeros(src.rows, src.cols);
		recon_src = decoded_error + epi_recon;
		double max_intensity = .0, min_intensity = .0;
		cv::minMaxLoc(recon_src, &min_intensity, &max_intensity);
		if ((min_intensity < 0) | (max_intensity > 255)){/*If Intensities overflow*/
			std::cout << "Intensities overflow!\n";
			// intensity should belong to the range [0, 255];
			for (int i = 0, rowNum = recon_src.rows; i != rowNum; ++i)
				for (int j = 0, colNum = recon_src.cols; j != colNum; ++j){
				recon_src.at<double>(i, j) = std::max<double>(recon_src.at<double>(i, j), 0.0);// >= 0;
				recon_src.at<double>(i, j) = std::min<double>(recon_src.at<double>(i, j), 255.0);// <= 255
				}
	}/*End of If Intensities overflow*/

#ifdef _DEBUG
		if (!epi_recon.data)
			cout << "Mat epi_recon is NULL!" << endl;
#endif
#ifdef _DEBUG
		if (recon_src.data){
			cout << "Mat recon_src (i.e., Mat_<double> recon_src = decoded_error + epi_recon) is:\n";
			display_n_by_n_elements(recon_src, 10);
		}
		else
			cout << "Mat recon_src (i.e., Mat_<double> recon_src = decoded_error + epi_recon) is NULL!\n";
#endif
		if (IsSaveImgFlag){
		bool IsWtirren = imwrite(recon_img_name, recon_src);
#ifdef _DEBUG
		if (!IsWtirren)
			cout << "Cannot write the image : " << recon_img_name << endl;
#endif
		}

		getEpiPSNR_MSE(src, recon_src, psnr, mse);
		if (fs->is_open()){
			*fs << "Uniform Quantization level = " << s_jpeg_quality << ", nlevels = " << uniQuant_BINARY_FILE.nLevels
				<< ", bits = " << bit_for_quantized_error << ", psnr = " << psnr << ", mse = " << mse << endl;
#ifdef _DEBUG
			std::cout << "Operation successfully performed\n";
#endif
			
		}
		else{
			std::cout << "Error opening std::fstream file";
		}
	} /*end of decoding*/

	delete [] c_jpeg_quality;
	delete [] char_error_binary_file_path;

}


// overload + 1
void Distortion_Bpp::doBMPUniformQuantDecoder(
	// const int & jpeg_quality, //  quantization level belonging to [1, 100], and higher values mean higher resolution;
	const Mat_<double> & src, // original image
	const Mat_<double> & epi_recon, // epitome reconstruction
	const string & binary_error_save_path, // indicating binary file for saving the quantization result
	const Mat_<uchar> & quantized_error, //  Uniform quantization result of residual images in the BMP format.
	// here we use uchar, since we will guarantee that the quantized value level is equal to or less than 256.
	const string & recon_img_name, // final reconstruction result from epitome reconstruction and uniform quantization of residual images.
	const bool & IsSaveImgFlag,
	double & psnr,
	double & mse,
	std::fstream * fs
	){

	//********************************************
	// some variables for uniform quantization
	//********************************************
	char * c_jpeg_quality = new char[10];
	// char c_jpeg_quality[10];
	// set three-character length
	sprintf(c_jpeg_quality, "%03d", jpeg_quality);
	string s_jpeg_quality(c_jpeg_quality);

	// Function : FILE * fopen ( const char * filename, const char * mode );
	FILE * p_file = NULL;
	BitStream bs(p_file);
	UniformQuant uniQuant_BINARY_FILE;

	
	// syntax: char * strcpy(char * destination, const char * source);
	// Copies the C string pointed by source into the array pointed by destination, 
	// including the terminating null character (and stopping at that point).
	char * char_error_binary_file_path = new char[binary_error_save_path.length() + 1];
	std::strcpy(char_error_binary_file_path, binary_error_save_path.c_str());

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
		int bit_for_quantized_error;
		if (uniQuant_BINARY_FILE.nLevels == 1)
			bit_for_quantized_error = 1;
		else
			bit_for_quantized_error = (int)ceil(log2(uniQuant_BINARY_FILE.nLevels));

		//************************************
		// for body
		//************************************
		Mat_<double> decoded_error(src.rows, src.cols);
		uniQuant_BINARY_FILE.deUniformQuantize(quantized_error, decoded_error);
		
		Mat_<double> recon_src(src.rows, src.cols);
		recon_src = decoded_error + epi_recon;
		double max_intensity = .0, min_intensity = .0;
		cv::minMaxLoc(recon_src, &min_intensity, &max_intensity);
		if ((min_intensity < 0) | (max_intensity > 255)){/*If Intensities overflow*/
			std::cout << "Intensities overflow!\n";
			// intensity should belong to the range [0, 255];
			for (int i = 0, rowNum = recon_src.rows; i != rowNum; ++i)
				for (int j = 0, colNum = recon_src.cols; j != colNum; ++j){
				recon_src.at<double>(i, j) = std::max<double>(recon_src.at<double>(i, j), 0.0);// >= 0;
				recon_src.at<double>(i, j) = std::min<double>(recon_src.at<double>(i, j), 255.0);// <= 255
				}
		}/*End of If Intensities overflow*/

		if (IsSaveImgFlag){
		bool IsWtirren = imwrite(recon_img_name, recon_src);
#ifdef _DEBUG
		if (!IsWtirren)
			cout << "Cannot write the image : " << recon_img_name << endl;
#endif
		}
		getEpiPSNR_MSE(src, recon_src, psnr, mse);
		if (fs->is_open()){
			*fs << "Uniform Quantization level = " << s_jpeg_quality << ", nlevels = " << uniQuant_BINARY_FILE.nLevels
				<< ", bits = " << bit_for_quantized_error << ", psnr = " << psnr << ", mse = " << mse << endl;
#ifdef _DEBUG
			std::cout << "Operation successfully performed.\n";
#endif		
		}
		else{
			std::cout << "Error opening std::fstream file";
		}
	} /*end of decoding*/

	delete[] c_jpeg_quality;
	delete [] char_error_binary_file_path;

}