#ifndef _LearnEpitome_HPP
#define _LearnEpitome_HPP
#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <fstream>
#include <string>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <complex> // C++, not C
#include <fftw3.h>
#include <math.h>
#include <algorithm> 
#include <cstdlib>     /* srand, rand */
#include <ctime>       /* time */
#include <cfloat> // DBL_MAX = 1E+37 or greater
#include <tuple> // e.g., std::tuple, std::get, std::tie, std::ignore
#include <vector>

#include <random>  // a random number generator
#include "makeDirectory.h" // read images from some directory
#include "readWriteDataToFromFile.hpp" 
#include "filesStorage.hpp"

// FFTW_Convolution or Correlation
#include "convolution_fftw.h"

typedef double Epitome_Doulbe;
typedef float Epitome_Float;
#define tolerance 1e-10 // if less than this tolerance, it will be recognized as ZERO.
#define minVar 1e-6 // enforce a minimum variance.
// #define img_read_flag 1
// maximum in double

#define DOUBLE_MAX  DBL_MAX // no semicolons needed, since it is not a statement but a directive.
// minimum in double
#define DOUBLE_MIN  -1*DOUBLE_MAX // no semicolons needed, since it is not a statement but a directive.

using namespace std;
using namespace cv;
namespace c_fft = FFTW_Convolution;

typedef enum
{
	ORIGINAL_FFTW,
	UPDATED_FFTW,
	NAME_DIFFERENCE,
	IMG_ENCODE_TYPE,
	OVERHEAD_SIZE
} EpitomeVariableName;

// ImgEpiotme is able to learn epitome based on both color images and gray-scale images.
class ImgEpitome{
private:
	int NUM_THREADS; // thread numbers used for the parallel calculation;
	int MAX_PARALLEL_SIZE; // due to the memory limit, we do parallel computation several times,
	// for each parallel computation, the OMP for(int i =0; i < MAX_PARALLEL_SIZE; ++ i)-loop 
	// can accept MAX_PARALLEL_SIZE as its ternimation limit;
	const int eWidth; //  # of columns
	const int eHeight; //# of rows
	const int patchSideLengh; // size of patches, with the width same as the height.
	int patchSpacing; // the distance between sampled patches from the input image.  it should be >= 2
	const bool isgrayScale; // true if gray images, or false if color images;
	const bool verbose; // for display the process of EM and each training patch by calculating the approximate left time.
	const bool RandomShifting; // random wiggle of the fixed xPatchIdx, this parameter is set to provide convenience to debugging
	unsigned int  Read_eMean_via_File_Flag; // read the v_eMean from ".yml"f files or images, like ".png" or ".jpg" images;
	const int numIteration;// for Expectation-Maximization Algorithm Iteration
	const std::string DatabaseDir;
	const std::string EpitomeResultDir;
	const std::string ReconsCompresDir;
	const std::string whatkindImgs;
	int Epi_Image_Read_Flag;
	std::string nameDifference;// e.g., = = "recon_"
	std::string imgEncodeType; // e.g. = ".jpg", ".bmp", or ".png";
	int overhead; // file size of epitome overhead, Units = Bytes

	//Epitome_Doulbe * v_eMean; //the Epitome, eMean, all zeros as initialized values
	//Epitome_Doulbe * v_eVar; // the Epitome, eVar , all ones as initialized values
	vector<Epitome_Doulbe> v_eMean_Red; //the Epitome, eMean of Red channel, all zeros as initialized values
	vector<Epitome_Doulbe>  v_eVar_Red; //the Epitome, eVar of Red channel, all ones as initialized values
	
	vector<Epitome_Doulbe> v_eMean_Gre; //the Epitome, eMean of Green channel, all zeros as initialized values
	vector<Epitome_Doulbe>  v_eVar_Gre; //the Epitome, eVar of Green channel , all ones as initialized values
	
	vector<Epitome_Doulbe> v_eMean_Blu; //the Epitome, eMean of blue channel, all zeros as initialized values
	vector<Epitome_Doulbe>  v_eVar_Blu; //the Epitome, eVar of blue channel, all ones as initialized values
	
	// for gray-scale images
	// for convenience, we allocate the variables for gray-scale images,
	vector<Epitome_Doulbe> v_eMean_Gray; //the Epitome, eMean of Red channel, all zeros as initialized values
	vector<Epitome_Doulbe>  v_eVar_Gray; //the Epitome, eVar of Red channel, all ones as initialized values

	vector<vector<int>> v_maxPost_row_column_Idx; // max-posterior indices read from the beforehand saved files (e.g.".yml")
	
	//	vector<vector<int>> v_row_Idx; // row indices of xPatches
	//	vector<vector<int>> v_col_Idx; // column indices of xPatches


	//****************************************************************************
	//****************************************************************************
	//****************************************************************************
	// all the variables for the parallel computation    *************************
	// including the shared ones, ************************************************
	// and the container for saving all the private results from each thread; ****
	//****************************************************************************
	// all the variables, which have nothing to do with xPatch (R, G, B or Gray),
	// should be considered as the shared one, i.e., 
	// the global variables which can be used by all the parallel threads;
	
	// shared variables for the each xPatch (R, G, B or Gray) 
	// used in the parallel calculation;
	// Src Signal: epitome(R,G,B channels, or gray-scale)

	/*constants*/
	vector<Epitome_Doulbe> v_patchSize_Ones;
	/*R channel*/
	vector<Epitome_Doulbe> v_eMeanOverVar_R;
	vector<Epitome_Doulbe> v_InvVar_R;
	vector<Epitome_Doulbe> v_LogVar_R; // not xPatch
	vector<Epitome_Doulbe> v_eMean2OverVar_R; // not xPatch
	vector<Epitome_Doulbe> v_eLogVarSum_R;
	vector<Epitome_Doulbe> v_eeSum_R;
	/*G channel*/
	vector<Epitome_Doulbe> v_eMeanOverVar_G;
	vector<Epitome_Doulbe> v_InvVar_G;
	vector<Epitome_Doulbe> v_LogVar_G; // not xPatch
	vector<Epitome_Doulbe> v_eMean2OverVar_G; // not xPatch
	vector<Epitome_Doulbe> v_eLogVarSum_G;
	vector<Epitome_Doulbe> v_eeSum_G;
	/*B channel*/
	vector<Epitome_Doulbe> v_eMeanOverVar_B;
	vector<Epitome_Doulbe> v_InvVar_B;
	vector<Epitome_Doulbe> v_LogVar_B; // not xPatch
	vector<Epitome_Doulbe> v_eMean2OverVar_B; // not xPatch
	vector<Epitome_Doulbe> v_eLogVarSum_B;
	vector<Epitome_Doulbe> v_eeSum_B;
	/*Gray channel*/
	vector<Epitome_Doulbe> v_eMeanOverVar_Gray;
	vector<Epitome_Doulbe> v_InvVar_Gray;
	vector<Epitome_Doulbe> v_LogVar_Gray; // not xPatch
	vector<Epitome_Doulbe> v_eMean2OverVar_Gray; // not xPatch
	vector<Epitome_Doulbe> v_eLogVarSum_Gray;
	vector<Epitome_Doulbe> v_eeSum_Gray;

	// those values are the final results,
	// since we do parallel calculation for all the patches extracted from one image;
	// e.g., for image 1, thus the temporary summation of sumQ, sumQX, or sumQXX is
	// just from this image; but the real sumQ, sumQX (R, G, B, or gray) and 
	// sumQXX (R, G, B, or gray), are gained through the summation of all the images;

	// before the starting of extracted training patches loop, 
	// please make sure zero-initialized doubles to them.
	// i.e., clear the matrices used for collecting sufficient statistics.
	// zero-initialized doubles
	vector<Epitome_Doulbe> v_double_sumQ,
		v_double_sumQX_R, v_double_sumQXX_R,
		v_double_sumQX_G, v_double_sumQXX_G,
		v_double_sumQX_B, v_double_sumQXX_B,
		v_double_sumQX_Gray, v_double_sumQXX_Gray;

	// for saving the max_unnormalized_post for each image patch out of one image;
	vector<Epitome_Doulbe> v_max_unnormalized_post;

	/*
	// all the variables, which involve with xPatch (R, G, B or Gray),
	// should be considered as private ones, i.e., 
	// the private variables which can be used by some thread (i.e., one thread);
	//R channel
	// those variables should have (PatchNum = r_size* c_size) elements;
	vector<Epitome_Doulbe * > v_p_sumQX_R;
	vector<Epitome_Doulbe * > v_p_sumQXX_R;
	//G channel
	// those variables should have (PatchNum = r_size* c_size) elements;
	vector<Epitome_Doulbe * > v_p_sumQX_G;
	vector<Epitome_Doulbe * > v_p_sumQXX_G;
	//B channel
	// those variables should have (PatchNum = r_size* c_size) elements;
	vector<Epitome_Doulbe * > v_p_sumQX_B;
	vector<Epitome_Doulbe * > v_p_sumQXX_B;

	//Gray channel
	// those variables should have (PatchNum = r_size* c_size) elements;
	vector<Epitome_Doulbe * > v_p_sumQX_Gray;
	vector<Epitome_Doulbe * > v_p_sumQXX_Gray;
	//Some variables is the summation of all the channels
	vector<Epitome_Doulbe * >  v_p_normalized_post; // normalized posterior
	vector<Epitome_Doulbe * >  v_p_sumQ;
	*/
	//****************************************************************************
	//****************************************************************************
	//****************************************************************************
	//****************************************************************************

public:
	// constructor
	ImgEpitome(int width, int height, int numImgs, int patchLengh, int patchSpac, int numItera, 
		bool isGray,  bool display, bool random,
		std::string DatabaseDir,
		std::string EpitomeResultDir,
		std::string ReconsCompresDir,
		std::string whatkindImgs,
		std::string nameDifference,// e.g., = = "recon_"
		std::string imgEncodeType, // e.g. = ".jpg", ".bmp", or ".png";
		int overhead, // file size of epitome overhead, Units = Bytes
		int nthreads, // thread numbers used for the parallel calculation;
		int max_omp_for_idx
		);  
	
	/*void randomizeEpitomeMean(
		Mat & mat_currentInputImage, // the current input image
		std::vector<double> & v_eMean //the Epitome, eMean, all zeros as initialized values
		);*/

	void normalize(Epitome_Doulbe * input, 
		int num, double minVal, double maxVal, double newMedian, double newMaxVal);
	
	void initialImgEpitome(vector<string> & v_ImgPath, // just image names (including file extension), without full directory in this parameter
		double weight_for_sigma, // e.g., = 1/100 
		// (see "plus Gaussian noise with 1/100th of the standard deviation (i.e., sigma) of the training data")
		bool IsConstantValueInitial,
		double constantVal
		);

	void randomEpitomeMean(vector<string> & v_ImgPath, // image names, including full directory in this parameter
		Epitome_Doulbe weight_for_sigma // e.g., = 1/100 
		// (see "plus Gaussian noise with 1/100th of the standard deviation (i.e., sigma) of the training data")
		);

	void updateEpitome(const char * value, EpitomeVariableName type);
	void clearImgEpitome();
	void displayEpitome();
	void parseString(const string & s_input, // e.g., == "4-8-12"
		const char & flag, // e.g., = '-'
		int * i_output // output which have been parsed based on the input string
		);

	// see http://stackoverflow.com/questions/9370493/inline-function-members-inside-a-class
	inline int ImgEpitome::getEpitomeWidth(){ return eWidth; }
	inline int ImgEpitome::getEpitomeHeight(){ return eHeight; }
	inline int ImgEpitome::getEpitomePatchSideLength(){ return patchSideLengh; }
	inline int ImgEpitome::getEpitomePatchSpacing(){ return patchSpacing; }
	inline int ImgEpitome::getEpitomeNumIteration(){ return numIteration; }
	inline bool ImgEpitome::getEpitomeVerbose(){ return verbose; }
	inline bool ImgEpitome::getEpitomeRandomShifting(){ return RandomShifting; }
	inline bool ImgEpitome::getEpitomeIsGrayScale(){ return isgrayScale; }
	inline void ImgEpitome::setFileExtension(std::string & filename){
		// make sure the filename is imencodeType, say, ".jpg"
		std::size_t pos = filename.find(".");      // e.g., position of ".jpg" in nameTemp
		filename = filename.substr(0, pos) + imgEncodeType;     // get beginning to ".jpg"
	}

	// rand() % 2 : generate a random integer within [0, 1], due to the operation "%".
	// int rand(void) : Returns a pseudo-random integral number in the range between 0 and RAND_MAX.
	// (double)rand() / RAND_MAX : generate a random double floating variable within [0, 1]
	// NOTE: how to generate random numbers between two doubles in c++ 
	// double fRand(double fMin, double fMax)
	// {
	//	double f = (double)rand() / RAND_MAX;
	//	return fMin + f * (fMax - fMin);
	// }
	// And remember to call srand() with a proper seed each time your program starts.

	//	v_eMean.at(temp) = ((double)rand() / RAND_MAX) * pixelStd + pixelMean;
	//	v_eMean.at(temp) = nd(normal_seed); 
	//	v_eMean.at(temp) = (rand() % 2); //*pixelStd1 + pixelMean; // rand()% 2 can only generate random integers, due to the operation "%".
	inline void getNormal_distribution( 
		vector<Epitome_Doulbe> & v,
		double & mu, // mean
		double & sigma, // standard deviation
		double & minCoeff,
		double & maxCoeff,
		double & weight_for_sigma // e.g., = 1/100 
		// (see "plus Gaussian noise with 1/100th of the standard deviation (i.e., sigma) of the training data")
		){
		// to initialize the v_eMean using the input image pixel mean and standard deviation.
		srand(time(NULL)); //Initialize random number generator
		std::default_random_engine normal_seed(time(NULL)); //seed
		std::normal_distribution<double> nd(mu, weight_for_sigma *sigma); //mean followed by a standard deviation
		
		// pay attention the initialization value here
		minCoeff= DOUBLE_MAX;
		maxCoeff = DOUBLE_MIN;
		for (int r = 0; r < eHeight; ++r){ /*row*/
			for (int c = 0; c < eWidth; ++c){ /*column*/
				int temp = r*eWidth + c;
				v[temp] = nd(normal_seed);
				minCoeff = minCoeff < v[temp] ? minCoeff : v[temp];
				maxCoeff = maxCoeff > v[temp] ? maxCoeff : v[temp];
			}/*end of row*/
		} /*end of column*/
		
		std::cout << "minCoeff =  " << minCoeff << " , and maxCoeff =  " << maxCoeff << endl;
	}


	tuple<vector<int>, vector<int>> getRowColIdxOfOneImg(
		const int & image_h, const int & image_w);

	void ImgEpitome::learnEpitome(EpitomeVariableName flag);
	
	void ImgEpitome::reconImgsAfterLearning(
		const int & patchSpace, 
		const unsigned int & Read_eMean_Flag);

	// this function is  made to reconstruct images based on the variant row-column indices
	// which are derived from the baseline epitome result via quantization and/or down- or up-sampling
	// to the baseline row-column-indices;
	void ImgEpitome::reconImgsDirect(
		const int & patchSpace,
		const string & baselineEpitomeResultDir, // baseline Epitome result directory
		const unsigned int & Read_eMean_Flag);

	// FFT calculation via FFTW, version 1
	void ImgEpitome::learnImgGenreEpitome(
		const vector<string> & v_imageFileNameList // read input images via their file-names.
		);
	
	// FFT calculation via FFTW, version 2,
	// specifically, we use the function introduced by the file "convolution_fftw.h" 
	// and "convolution_fftw.cpp";
	void ImgEpitome::learnImgGenreEpitome_plus(const vector<string> & v_imageFileNameList);

	void ImgEpitome::getMaxPostRowColIdx(
		const vector<string> & v_imageFileNameList, // read input images via their file-names.
		const int & patchSpac
		);


	// build reconstructed images based on the learned Epitome and indices, without the input images
	void ImgGenreReconViaIdx(const vector<string> & v_imageFileNameList // the names of reconstructed images.
		);

	bool ImgEpitome::get_xPatch_unnormalPosterior_RGB_or_Gray_Channel(
		// It's pretty much all written in the FFTW documentation about thread safety:
		// but some care must be taken because the planner routines share data
		// (e.g.wisdom and trigonometric tables) between calls and plans.
		// so here we set ws as a parameter,
		// making those threads share the some fftw plan in ws,
		// in a typical application of FFT, we create one plan, 
		// and execute it several times, since the function fftw_execute (...)
		// is the only thread-safe (re-entrant) routine in FFTW.
		c_fft::FFT_Workspace & ws,
		// used saving maxPost_row_column_Idx, for identifying each patch idx;
		// to save the maxPost_row_column_Idx or not, 
		// since we just want to save the indices during the last time EM iteration;
		const int  patchIdx,
		const bool IsSaveRowColIdx,
		const int  temp_location, //  temp_location = temp_r * image_w + temp_c;
		const int  image_w,
		const int  current_image_index, // used for saving maxPost_row_column_Idx;
		const cv::Mat & mat_currentInputImage,
		const int  PatchNum, // used for disply the remaining time;
		Epitome_Doulbe * const p_max_unnormalized_post
		);


	bool ImgEpitome::get_xPatch_sumQXX_etc_RGB_or_Gray_Channel(
		// It's pretty much all written in the FFTW documentation about thread safety:
		// but some care must be taken because the planner routines share data
		// (e.g.wisdom and trigonometric tables) between calls and plans.
		// so here we set ws as a parameter,
		// making those threads share the some fftw plan in ws,
		// in a typical application of FFT, we create one plan, 
		// and execute it several times, since the function fftw_execute (...)
		// is the only thread-safe (re-entrant) routine in FFTW.
		c_fft::FFT_Workspace & ws,
		const int  patchIdx, // used saving maxPost_row_column_Idx, for identifying each patch idx;
		// to save the maxPost_row_column_Idx or not, 
		// since we just want to save the indices during the last time EM iteration;
		const bool IsSaveRowColIdx,
		const int  temp_location, //  temp_location = temp_r * image_w + temp_c;
		const int  image_w,
		const int  current_image_index, // used for saving maxPost_row_column_Idx;
		const cv::Mat & mat_currentInputImage,
		const int  PatchNum, // used for disply the remaining time;
		Epitome_Doulbe * const p_max_unnormalized_post,
		// v_p_sumQ_sumQX_RGB_etc[7]  = { 
		//  p_sumQ, // summation of 3 channels;
		//  p_sumQX_R, p_sumQXX_R, 
		//  p_sumQX_G, p_sumQXX_G,
		//  p_sumQX_B, p_sumQXX_B,
		//}; 
		vector<Epitome_Doulbe * const> v_p_sumQ_sumQX_RGB_etc
		);

	void ImgEpitome::paral_getMaxPostRowColIdx(
		const vector<string> & v_imageFileNameList, // read input images via their file-names.
		const int & patchSpac);

	vector<vector<int>> read_row_col_idx_form_YML(
		const string & s_Idx //  e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
		);

	void write_row_col_idx_2_YML(
		const string & s_dst_Idx, //  e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
		int & numInputImage // the number of images
		);
		
	void encodeMaxRowColIdx();
	void decodeMaxRowColIdx();

};
#endif