#ifndef _ROW_COL_IDX_TO_J2K_IMAGE_H
#define _ROW_COL_IDX_TO_J2K_IMAGE_H
// change row_column indices into j2k image for a less compressed file size
#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include <vector>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "iostream"
#include "fstream"
#include "makeDirectory.h"
#include "quantizer.h"
#include "ShowManyImages.h"

// boost filesystem
#define BOOST_FILESYSTEM_NO_DEPRECATED 
// For new code, defining "BOOST_FILESYSTEM_NO_DEPRECATED" before including filesystem headers is strongly recommended. 
// This prevents inadvertent use of old features, particularly legacy function names, 
// that have been replaced and are going to go away in the future.
#include <boost/filesystem/operations.hpp>
using namespace boost::filesystem;
using namespace std;
using namespace cv;
namespace bf = boost::filesystem;
#define MAX_CHAR_NUM_OF_FILES_PATH  800
#define MAX_VALUE_FOR_UCHAR 256
#define INTERPOLATION_METHOD cv::INTER_CUBIC
/*
INTER_NEAREST - a nearest-neighbor interpolation
INTER_LINEAR - a bilinear interpolation (used by default)
INTER_AREA - resampling using pixel area relation. It may be a preferred method for image decimation, 
             as it gives moire?free results. But when the image is zoomed, it is similar to the INTER_NEAREST method.
INTER_CUBIC - a bicubic interpolation over 4x4 pixel neighborhood
INTER_LANCZOS4 - a Lanczos interpolation over 8x8 pixel neighborhood
*/


vector<vector<int>> read_row_col_idx_form_YML(
	const string & whatkindImgs, // e.g., = "Edwin-less"
	const string & s_Idx, //  e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
	int & numInputImage
	);

void write_row_col_idx_2_YML(
	const string & whatkindImgs, // e.g., = "Edwin-less"
	const string & s_dst_Idx, //  e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
	const vector<vector<int>> v_dst_maxPost_row_column_Idx, // indices to be written into the YML files
	int & numInputImage // the number of images
	);

void getRowColIdxSizeOfOneImg(const int & patchSideLengh, // = 8
	const int & patchSpacing, // = 4, or 8
	const int & image_h, // image height
	const int & image_w, // image width
	const bool & verbose,
	int & r_idx_size,
	int & c_idx_size,
	vector<int> & v_rIdx, 
	vector<int> & v_cIdx
	);

//*********************************************************************************************************
// this function is used to get the size (unit: Bytes) of all the files (like, images) in some a directory
//*********************************************************************************************************
uintmax_t getFileSizeFromDir(const string & input_img_dir // input image category directory
	);

// 1/4 down-sample of mat  
// recon = interpolated  up-sample of the 1/4 down-sample result,
// that is , the dimension of recon is the same as that of src.
cv::Mat quarter_down_sample_mat(cv::Mat & src, cv::Mat & recon);

cv::Mat col_half_down_sample_mat(cv::Mat & src, cv::Mat & recon);

// change row_column indices into j2k image for a less compressed file size
uintmax_t saveRowColIdx2_J2KImg(
	const int & patchSideLength,
	const int & patchSpacing,
	const int & eWidth,
	const string & whatkindImgs,
	const string & JPEG2K_EXE_Base_Path,
	const string & category_RowColIdx, // e.g., "Edwin-less_4_0_RowCol_Idx";
	const string & path_Idx, // e.g., = ".../Edwin-less_4_RowCol_Idx.yml"
	const string & J2K_Result_BasePath, //
	const bool & is7zCompression, // do or not do compression from bmp files to 7z files;
	uintmax_t & size_7z_file // if do 7z compression, then find its file size; otherwise set the value of 0;
	);

void down_up_sampling_2_row_col_Idx(
	const int & r_idx_size,
	const int & c_idx_size,
	const int & imgIdx,
	const int & image_h,
	const int & image_w,
	const int & eWidth,
	const bool & IsShowingImgFlag,
	const vector<int> & v_maxPost_row_Idx,
	const vector<int> & v_maxPost_col_Idx,
	const double & f_width, // scale factor along the horizontal axis; it is computed as (double)dsize.width / src.cols;
	const double & f_height,// scale factor along the vertical axis; it is computed as (double)dsize.height / src.rows;
	vector<vector<int>> & v_dst_maxPost_row_column_Idx
	);

void quan_dequantize_row_col_Idx(
	const int & r_idx_size,
	const int & c_idx_size,
	const int & eWidth,
	const double QuantizationLevel, // constant quantization level value, 0 <= QuantizationLevel <= 100;
	vector <int> & v_maxPost_row_Idx,
	vector<int> & v_maxPost_col_Idx,
	vector <int> & v_quantized_maxPost_row_Idx, // quantized row indices;
	vector <int> & v_quantized_maxPost_col_Idx // quantized column indices;
	);

// do uniform quantization to row-column indices of epitome
// input : existing indices in the format of YML files;
// output: quantized or even downsampled indices, still in the format of YML;
void UniformQuant_row_col_idx_saved_as_YML(
	const string & whatkindImgs, // e.g., = "Edwin-less"
	const string & s_src_Idx, //  path of existing input indices, e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
	const string & s_dst_Idx, //  path of output indices, e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx_plus.yml";
	const string & s_quantized_Idx,
	const int & patchSideLength,
	const int & patchSpacing,
	const int & eWidth,
	const bool & IsShowingImgFlag,
	const bool & IsUniformQuantize, // uniform quantization flag;
	const double QuantizationLevel, // constant quantization level value, 0 <= QuantizationLevel <= 100;
	const bool & IsDownSampled, // down-sampling or up-sampling flag;
	const double & f_width, // scale factor along the horizontal axis; it is computed as (double)dsize.width / src.cols;
	const double & f_height// scale factor along the vertical axis; it is computed as (double)dsize.height / src.rows;
	);



#endif