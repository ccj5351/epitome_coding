#ifndef _DISTORTION_BIT_PER_PIXEL_H
#define _DISTORTION_BIT_PER_PIXEL_H

// see http://stackoverflow.com/questions/14448433/error-c4996-fopen-not-declared
#define _CRT_SECURE_NO_DEPRECATE
#include <cstdio> // C header

#include <iostream>
#include <map>
#include <fstream>
#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp> // histogram
#include <opencv2/highgui/highgui.hpp>
#include <math.h>
#include <algorithm> 
#include <cstdint> // uintmax_t
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <cfloat> // DBL_MAX = 1E+37 or greater
#include <tuple> // e.g., std::tuple, std::get, std::tie, std::ignore
#include <vector>
#include <random>  // a random number generator
#include <stdio.h>
#include <stdlib.h>


// boost filesystem
#define BOOST_FILESYSTEM_NO_DEPRECATED 
// For new code, defining "BOOST_FILESYSTEM_NO_DEPRECATED" before including filesystem headers is strongly recommended. 
// This prevents inadvertent use of old features, particularly legacy function names, 
// that have been replaced and are going to go away in the future.
#include <boost/filesystem/operations.hpp>
using namespace boost::filesystem;
namespace bf = boost::filesystem;

// /* There is a better solution, if you have the Boost C++ library installed which provides us with the following function:
//  * boost::uintmax_t file_size(const path& p); 
//  */


#include "makeDirectory.h" // read images from some directory 
#include "bitstream.h" // bit stream writing to and reading from physical files;
#include "quantizer.h" // uniform quantization 

// maximum in double
#define DOUBLE_MAX  DBL_MAX // no semicolons needed, since it is not a statement but a directive.
// minimum in double
#define DOUBLE_MIN  -1*DOUBLE_MAX // no semicolons needed, since it is not a statement but a directive. 

// #define Epi_Image_Read_Flag 1
#define MAX_CHAR_NUM_OF_FILES_PATH 800
#define MAX_DATE 80 /* for strftime(..) function!*/
// #define SAVE_IMG_NUM 10 /*the number of images will be saved as physical files*/

using namespace std;
using namespace cv;
#define  BYTE  unsigned char
#define uint unsigned int

typedef enum
{
	JPEG,
	JPEG2000,
	PNG,
	BMP
} EncodeTypeName;

typedef enum
{
	SET2MEDIAN,
	BETWEEN_MIN_MAX,
	NOTHING,
	BIT_PLANE_SLICING_COMBINE4BITSIMG_JPEG,
	BIT_PLANE_SLICING_JPEG,
	UNIFORM_QUANTIZATION_JPEG,
	UNIFORM_QUANTIZATION_BINARY_FILE
} QuantizeType;

typedef enum {
	UNIFORM_QUANTIZATION, // 0
	UNIFORM_QUANTIZATION_LOSSLESS_J2K, // 1
	UNIFORM_QUANTIZATION_LOSSY_JPEG, // 2
	UNIFORM_QUANTIZATION_LOSSY_J2K, // 3
} COMPRE_FLAG_2_EPI_RECON;


typedef enum {
	plane_1,
	plane_2,
	plane_3,
	plane_4,
	plane_5,
	plane_6,
	plane_7,
	plane_8
} e_bit_plane; // e_ : enum

typedef struct{
	BYTE mask;
	string s_bit_plane;
	uint i_th_bit;
} s_bit_plane; // s_ ; struct

// implement standard codecs like, JPEG, JPEG2000, etc compressor to epitome reconstructed images
class Distortion_Bpp{
private:
	string DatabaseDir ;
	string EpitomeResultDir;
	string ReconsCompresDir;
	EncodeTypeName encodeType; // e.g., jpeg, jpeg2000, bmp, png, etc.
	const string s_encodeType;
	QuantizeType quantizeType; 
	string s_quantizeType;
	const int errorImgBit ;
	const int jpeg_quality; // cannot exceed the range of [0, 100];
	const string encodeTypeforReconsImgs; // e.g. = ".bmp";
	const string encodeTypeforErrorImgs; // e.g. = ".jpg";
	const string nameDifference; // e.g., = "recon-";
	const string whatKindofImgs;
	const int errHistThres;
	const int save_img_number;
	const bool isGrayScale; // true if gray images, or false if color images; 
	int img_read_flag;

public:
	// constructor
	Distortion_Bpp( string & input, string & epi, string & com, string & s_standard_encodeType, EncodeTypeName & encodeType, QuantizeType & quantizetype,
		string & s_quantizeType, int  & errImgBit, int & jpegQual,
		string & encodeTypeforReconsImgs,
		string & encodeTypeforErrorImgs,
		string & nameDifference,
		string & whatKindofImgs,
		int & errHistThres,
		int & save_img_num,
		bool & isGary
		);

	inline void setFileExtension(std::string & filename, const string & extension){
		// make sure the filename is imencodeType, say, ".jpg"
		std::size_t pos = filename.find(".");      // e.g., position of ".jpg" in nameTemp
		filename = filename.substr(0, pos) + extension;     // get beginning to ".jpg"
	}

	// to check if a file exist
	inline bool is_File_Exist(const std::string& name) {
		ifstream f(name.c_str());
		if (f.good()) {
			f.close();
			return true;
		}
		else {
			f.close();
			return false;
		}
	}
	
	/*
	// see http://stackoverflow.com/questions/758001/log2-not-found-in-my-math-h
	log2() is only defined in the C99 standard, not the C90 standard. 
	Microsoft Visual C++ is not fully C99 compliant (heck, there isn't a single fully C99 compliant compiler in existence, 
	I believe -- not even GCC fully supports it), so it's not required to provide log2().
	*/
	// Calculates log2 of number.
	inline  double log2( double n ){  
    // log(n)/log(2) is log2.  
    return log( n ) / log(2.0);  
}
	inline void setImgReadFlag(){
		if (isGrayScale)
			img_read_flag = 0;
		else
			img_read_flag = 1;
	}
	// to check if a file exist
	// overload + 1
	inline bool is_File_Exist(const char *fileName){
		std::ifstream f(fileName);
		if (f.good()) {
			f.close();
			return true;
		}
		else {
			f.close();
			return false;
		}
	}

	 template<typename T>
	 inline void display_n_by_n_elements(const Mat_<T> & m_in, const int & n){
		 for (int i = 0; i != n; ++i){
			 for (int j = 0; j != n; ++j){
				 cout << m_in.at<T>(i, j) << "  ";
			 }
			 cout<< "\n";
		 }
		 cout << "\n";
}

	 inline void display_n_by_n_uchar_elements(const Mat_<uchar> & m_in, const int & n){
		 for (int i = 0; i != n; ++i){
			 for (int j = 0; j != n; ++j){
				 printf("%03d ",m_in.at<uchar>(i, j));
			 }
			 cout << "\n";
		 }
		 cout << "\n";
	 }

	inline  double vector_average(const vector<double> &v)
	 {
		 vector<double>::size_type taille = v.size();
		 double sum = .0;
		 for (vector<double>::const_iterator it = v.begin(); it != v.end(); ++it)
			 sum += *it;
		 return sum / taille;
	 }

	//declaration of functions
	double Distortion_Bpp::getPSNR(const Mat& input, const Mat & recons);

	// this function is used to get the size (unit: Bytes) of all the files (like, images) in some a directory
	uintmax_t Distortion_Bpp::getFileSizeFromDir(const string & input_img_dir // input image category directory
		);


	void  Distortion_Bpp::getJPEGCompressToEpiRecons(bool & Ishistgram);

	void Distortion_Bpp::getImgHistogram(
		const Mat & src, // only one input image
		const Mat & dst, // output histogram of the only 1 input image
		const int & histSize, // the number of bins , e.g., histSize = 256
	//	const float ** ranges, // the ranges ( for B,G,R) 
		bool & uniform, //
		bool & CurveTypeHist //
		);

	tuple<int, int, int, double , double> Distortion_Bpp::getGrayImgHistogram(
		const Mat & src, // only one input image
		const int & histSize, // the number of bins , e.g., histSize = 256 
		const int & hist_w, // e.g., hist_w = 512 or 1024, width of histogram drawing
		const int & hist_h, // e.g., hist_h = 400; height of histogram drawing
		bool & IsCurveTypeHist, // Flag indicating to draw the histogram, in the form of curves or rectangles
		const string & pathHist, // save histogram
		const string & foutHistVal
		);


	void Distortion_Bpp::imgCompressViaStandardCodec(
		const string & GetJPEGCompressDir, // = "E:\\ComputerVisionCCJ\\Epitome\\imageDatabase\\simi-objects-bmpImages";
		const string & JPEG2K_EXE_Base_Path
		);

	void Distortion_Bpp::getStandardCodecRate_DistortionfromDir(
		const string & GetJPEGCompressDir // = "E:\\imageDatabase\\forPSNR\\compress-bmp-images\\compress-small-bmp-images";
		);

	tuple<double, double, double> Distortion_Bpp::getStandardCodecRate_Distortion(
		const string & GetJPEGCompressDir // reconstructed compressed image category directory
		);

	tuple<double, double, double, double> Distortion_Bpp::getEpiCompressRate_Distortion(
		const string & ReconsCompresDir,
		const uintmax_t & sum_input_size, // input image database file size (Units : Bytes)
		const uintmax_t & epitomeDataFileSize, // epitome Data FileSize, Units = Bytes
		const COMPRE_FLAG_2_EPI_RECON & compres_flag_2_epi_recon
		);

	void Distortion_Bpp::getEpiCompressRate_DistortionfromDir(
		std::map<string, uintmax_t> & m_epitomeDataFileSize, // epitome Data FileSize, Units = Bytes
		const uintmax_t & sum_input_size // input image database file size (Units : Bytes)
		);


	/*
	**************************************************************************
	*** Bit - plane slicing
	**************************************************************************
	*/

	cv::Mat Distortion_Bpp::get_i_th_BitPlane(const e_bit_plane & e_bp, s_bit_plane & s_bp, const cv::Mat & in_img);

	cv::Mat Distortion_Bpp::recon_via_4_MSB_BitPlane(const cv::Mat & in_img);
	cv::Mat Distortion_Bpp::recon_via_2_MSB_BitPlane(const cv::Mat & in_img);
	
	void Distortion_Bpp::getEpiPSNR_MSE(const Mat& input, const Mat& recons, double & psnr, double & rmse);

	void Distortion_Bpp::getEpiPSNR_MSE_FromDir(string & inputDir,
		string & epi_res_Dir,
		string nameDifference // = "recon-"
		);

	int Distortion_Bpp::UniformQuantize(const Mat_<double> & src, Mat_<unsigned int> & dst, UniformQuant & uniQuant);
	int Distortion_Bpp::deUniformQuantize(const Mat_<unsigned int> & src, Mat_<double> & dst, UniformQuant & uniQuant);

	int Distortion_Bpp::UniformQuantize(const Mat_<double> & src, Mat_<BYTE> & dst, UniformQuant & uniQuant);
	int Distortion_Bpp::deUniformQuantize(const Mat_<BYTE> & src, Mat_<double> & dst, UniformQuant & uniQuant);
	
	void Distortion_Bpp::doBinUniformQuantCodec(
		// const int & jpeg_quality, //  quantization level belonging to [1, 100], and higher values mean higher resolution;
		const Mat_<double> & src, // original image
		const Mat_<double> & epi_recon, // epitome reconstruction
		const bool isConstQuantizationLevel, // set a quantization level or not
		const int QuantizationLevel, // constant quantization level value
		const string & error_save_path, // indicating binary file for saving the quantization result
		// const unsigned int & encoding_decoding_flag, // encoding_decoding_flag = 0  means encoding; = 1 means decoding
		const string & recon_img_name, // final reconstruction result from epitome reconstruction and uniform quantization of residual images.
		const bool & IsSaveImgFlag,
		double & psnr,
		double & rmse,
		std::fstream * fs
		);

	void Distortion_Bpp::doBMPUniformQuantEncoder(
		// const int & jpeg_quality, //  quantization level belonging to [1, 100], and higher values mean higher resolution;
		const Mat_<double> & src, // original image
		const Mat_<double> & epi_recon, // epitome reconstruction
		const bool isConstQuantizationLevel, // set a quantization level or not
		const int QuantizationLevel, // constant quantization level value
		const string & binary_error_save_path, // indicating binary file for saving the quantization result
		const string & error_img_name //  Uniform quantization result of residual images in the BMP format.
		);

	void Distortion_Bpp::doBMPUniformQuantDecoder(
		const Mat_<double> & src, // original image
		const Mat_<double> & epi_recon, // epitome reconstruction
		const string & binary_error_save_path, // indicating binary file for saving the quantization result
		const string & error_img_name, //  Uniform quantization result of residual images in the BMP format.
		const string & recon_img_name, // final reconstruction result from epitome reconstruction and uniform quantization of residual images.
		const bool & IsSaveImgFlag,
		double & psnr,
		double & rmse,
		std::fstream * fs
		);

	// overload +1 
	void Distortion_Bpp::doBMPUniformQuantDecoder(
		const Mat_<double> & src, // original image
		const Mat_<double> & epi_recon, // epitome reconstruction
		const string & binary_error_save_path, // indicating binary file for saving the quantization result
		const Mat_<uchar> & quantized_error, //  Uniform quantization result of residual images in the BMP format.
		// here we use uchar, since we will guarantee that the quantized value level is equal to or less than 256.
		const string & recon_img_name, // final reconstruction result from epitome reconstruction and uniform quantization of residual images.
		const bool & IsSaveImgFlag,
		double & psnr,
		double & rmse,
		std::fstream * fs
		);

	void Distortion_Bpp::CompressReconToEpiRecons(
		const string JPEG2K_EXE_Base_Path, // JPEG2K compressor and de-compressor base path
		const COMPRE_FLAG_2_EPI_RECON compres_flag_2_epi_recon);

	void Distortion_Bpp::uniQuan_0(
		const cv::Mat & input,
		const cv::Mat & reconsInput,
		const string & binary_error_save_path,
		const string & final_recon_img_name,
		const bool & IsSaveImgFlag,
		double & psnr,
		double & mse,
		std::fstream * p_fs
		);
	void Distortion_Bpp::uniQuan_1(
		const cv::Mat & input,
		const cv::Mat & reconsInput,
		const string & JPEG2K_EXE_Base_Path,
		const string & bin_header_error_save_path,
		const string & j2k_error_save_path,
		const string & de_error_save_path,
		const string & error_save_path,
		const string & final_recon_img_name,
		const bool & IsSaveImgFlag,
		double & psnr,
		double & mse,
		std::fstream * p_fs
		);
	void Distortion_Bpp::uniQuan_2(
		const cv::Mat & input,
		const cv::Mat & reconsInput,
		const string & final_recon_img_name,
		const string & bin_header_error_save_path,
		const string & jpg_error_save_path,
		const string & error_save_path,
		const bool & IsSaveImgFlag,
		double & psnr,
		double & mse,
		std::fstream * p_fs
		);
	void Distortion_Bpp::uniQuan_3(
		const cv::Mat & input,
		const cv::Mat & reconsInput,
		const string & JPEG2K_EXE_Base_Path,
		const string & bin_header_error_save_path,
		const string & j2k_error_save_path,
		const string & de_error_save_path,
		const string & error_save_path,
		const string & final_recon_img_name,
		const bool & IsSaveImgFlag,
		double & psnr,
		double & mse,
		std::fstream * p_fs
		);
	void Distortion_Bpp::parseString(const string & s_input, // e.g., == "4-8-12"
		const char & flag, // e.g., = '-'
		int * i_output // output which have been parsed based on the input string
		);

	/*
	function : strftime, to format time as string;
	Syntax: size_t strftime (char* ptr, size_t maxsize, const char* format, const struct tm* timeptr );
	Parameters:
	  ptr : Pointer to the destination array where the resulting C string is copied.
	  maxsize : Maximum number of characters to be copied to ptr, including the terminating null-character.

	Copies into ptr the content of format, expanding its format specifiers into the corresponding values 
	that represent the time described in timeptr, with a limit of maxsize characters.
	*/

	inline std::string Distortion_Bpp::get_curr_date(void){

		time_t now;
		char the_date[MAX_DATE];
		the_date[0] = '\0';

		now = time(NULL);
		// %d	Day of the month, zero-padded (01-31), e.g., 23;
		// %m	Month as a decimal number(01 - 12), e.g., 08;
		// %Y	Year, e.g.,	2001;
		// %H	Hour in 24h format (00-23), e.g., 14;
		// %M	Minute (00-59), e.g., 55;
		// %S	Second (00-61), e.g., 02;

		if (now != -1){
			strftime(the_date, MAX_DATE, "%d_%H_%M_%S", gmtime(&now));
		}
		// one kind of string constructor -- string( char * );
		return std::string(the_date);
	}

};

#endif