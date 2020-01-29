#include "Idx2Img.h"
#include "ShowManyImages.h"


// encode the max-row-column indices
// DPCM -- Differential Pulse Code Modulation on DC components in JPEG
// DC component is large and varied, but often closed to previous value (like lossless JPEG)
// encode the difference from previous value

uintmax_t getFileSizeFromDir(const string & input_img_dir // input image category directory
	){
	uintmax_t sum_input_size = 0;
	vector<string> filelist;
	GetFileList(input_img_dir, &filelist);
	auto filelistSize = filelist.size();
	for (unsigned int i = 0; i < filelistSize; ++i){
		string input_path = input_img_dir + "\\" + filelist[i];
		sum_input_size += file_size(input_path);
	}
	vector<string>().swap(filelist);
	return sum_input_size;
}

vector<vector<int>> read_row_col_idx_form_YML(
	const string & whatkindImgs, // e.g., = "Edwin-less"
	const string & s_Idx, //  e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
	int & numInputImage
	){

	FileStorage fs_Idx_out(s_Idx, FileStorage::READ);
	vector<vector<int>> v_maxPost_row_column_Idx;
	if (!fs_Idx_out.isOpened()){
		cerr << "failed to open " << s_Idx << endl;
	}
	else{
		FileNode fn_Idx = fs_Idx_out[whatkindImgs];
		FileNodeIterator it = fn_Idx.begin(), it_end = fn_Idx.end();
		numInputImage = fn_Idx.size();
		int temp_idx = 0;
		v_maxPost_row_column_Idx = vector<vector<int>>(numInputImage, vector<int>());
		for (; it != it_end; ++it, temp_idx++){
			v_maxPost_row_column_Idx[temp_idx].clear(); // max-posterior indices read from the beforehand saved files (e.g.".yml")
			(*it)["maxPost_row_col_Idx"] >> v_maxPost_row_column_Idx[temp_idx];
		}
	}
	fs_Idx_out.release();
	return v_maxPost_row_column_Idx;

}

void write_row_col_idx_2_YML(
	const string & whatkindImgs, // e.g., = "Edwin-less"
	const string & s_dst_Idx, //  e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
	const vector<vector<int>> v_dst_maxPost_row_column_Idx, // indices to be written into the YML files
	int & numInputImage // the number of images
	){
	//******************************************************************
	// third, save the uniform quantized and/or down-sampled/up-sampled 
	// indices to the output YML file;
	//******************************************************************
	FileStorage fs_Idx_write(s_dst_Idx, FileStorage::WRITE);
	if (!fs_Idx_write.isOpened()){
		cerr << "failed to open " << s_dst_Idx << endl;
	}
	else{
		//*******************************************************
		// first write the beginning of the FileStorage file.
		//*******************************************************
		fs_Idx_write << whatkindImgs << "[:";

		/*each input image*/
		for (int current_image_index = 0; current_image_index != numInputImage; ++current_image_index){
			// save parameters into the ".YML" file for image reconstruction
			// 2-kind Parameters of current image include:
			//     * image size: image_h * image_w;
			//     * possible row and column indices of maximun posterior element of each xPatch in currently read image;
			fs_Idx_write << "{:"
				<< "maxPost_row_col_Idx" << v_dst_maxPost_row_column_Idx[current_image_index]
				<< "}";
		} /* end of each input image*/
	}
	// write the end of the FIleStorage file
	if (!fs_Idx_write.isOpened()){
		cerr << "failed to open " << s_dst_Idx << endl;
	}
	else
		fs_Idx_write << "]";
	fs_Idx_write.release();
}

// do quantization
void quan_dequantize_row_col_Idx(
	const int & r_idx_size,
	const int & c_idx_size,
	const int & eWidth,
	const double QuantizationLevel, // constant quantization level value, 0 <= QuantizationLevel <= 100;
	vector <int> & v_maxPost_row_Idx, // first to be quantized and then be de-quantized
	vector<int> & v_maxPost_col_Idx, // first to be quantized and then be de-quantized
	vector <int> & v_quantized_maxPost_row_Idx, // quantized row indices;
	vector <int> & v_quantized_maxPost_col_Idx // quantized column indices;
	){
	// do uniform quantization
	double max_row_col_Idx = eWidth - 1; // e.g., eWidth = 256, max_row_col_Idx = 255;
	double min_row_col_Idx = .0;
	int nLevels = (int)(QuantizationLevel * (max_row_col_Idx - min_row_col_Idx) / 100);
	// at last, make sure 1 <= nLevels <= 256
	nLevels = nLevels > MAX_VALUE_FOR_UCHAR ? MAX_VALUE_FOR_UCHAR : nLevels;
	nLevels = nLevels < 1 ? 1 : nLevels;
	UniformQuant uniQuant(max_row_col_Idx, min_row_col_Idx, nLevels);
	//vector <uint> v_quantized_maxPost_row_Idx(r_idx_size * c_idx_size, 0), v_quantized_maxPost_col_Idx(r_idx_size *c_idx_size, 0);
	// here we do not save the header of uniform quantization,
	uniQuant.UniformQuantize(v_maxPost_row_Idx, v_quantized_maxPost_row_Idx);
	uniQuant.UniformQuantize(v_maxPost_col_Idx, v_quantized_maxPost_col_Idx);
	// do de-quantization
	uniQuant.deUniformQuantize(v_quantized_maxPost_row_Idx, v_maxPost_row_Idx);
	uniQuant.deUniformQuantize(v_quantized_maxPost_col_Idx, v_maxPost_col_Idx);
	// release memory
	// vector<uint>().swap(v_quantized_maxPost_row_Idx);
	// vector<uint>().swap(v_quantized_maxPost_col_Idx);
}

// usage:
// do down-sampling and then up-sampling to the indices: v_maxPost_row_Idx, v_maxPost_col_Idx,
// and the result to the output of v_dst_maxPost_row_column_Idx;
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
	){
	
	// Size(double width, double height), parameters of down-sampling or up-sampling, i.e., output image size;
	cv::Size dst_size = cv::Size(c_idx_size*f_width,r_idx_size*f_height);
	Mat m_maxPost_row_Idx(r_idx_size, c_idx_size, CV_8UC1), m_maxPost_col_Idx(r_idx_size, c_idx_size, CV_8UC1),
		m_recon_row_Idx(r_idx_size, c_idx_size, CV_8UC1), m_recon_col_Idx(r_idx_size, c_idx_size, CV_8UC1),
		m_down_row_Idx(dst_size, CV_8UC1), m_down_col_Idx(dst_size, CV_8UC1);
	

	for (int i = 0; i < r_idx_size; ++i)
		for (int j = 0; j < c_idx_size; ++j)
			m_maxPost_col_Idx.at<uchar>(i, j) = v_maxPost_col_Idx[i * c_idx_size + j];

	for (int i = 0; i < r_idx_size; ++i)
		for (int j = 0; j < c_idx_size; ++j)
			m_maxPost_row_Idx.at<uchar>(i, j) = v_maxPost_row_Idx[i * c_idx_size + j];
	// firstly, do down-sampling
	// OpenCV class InputArray:
	// This is the proxy class for passing read-only input arrays into OpenCV functions. It is defined as
	// typedef const _InputArray& InputArray;
	// where _InputArray is a class that can be constructed from Mat, Mat_<T>, Matx<T, m, n>, std::vector<T>, std::vector<std::vector<T> > 
	// or std::vector<Mat>. It can also be constructed from a matrix expression.
	cv::resize(m_maxPost_row_Idx, m_down_row_Idx, cv::Size(), f_width, f_height, INTERPOLATION_METHOD);
	cv::resize(m_maxPost_col_Idx, m_down_col_Idx, cv::Size(), f_width, f_height, INTERPOLATION_METHOD);
	// secondly, do up-sampling for reconstruction
	cv::resize(m_down_row_Idx, m_recon_row_Idx, m_recon_row_Idx.size(), 0, 0, INTERPOLATION_METHOD);
	cv::resize(m_down_col_Idx, m_recon_col_Idx, m_recon_col_Idx.size(), 0, 0, INTERPOLATION_METHOD);

	if (IsShowingImgFlag){

		// Vector of Images
		std::vector<cv::Mat> v_Mat(6, m_maxPost_row_Idx);
		v_Mat[0] = m_maxPost_row_Idx;
		v_Mat[1] = m_maxPost_col_Idx;
		v_Mat[2] = m_down_row_Idx;
		v_Mat[3] = m_down_col_Idx;
		v_Mat[4] = m_recon_row_Idx;
		v_Mat[5] = m_recon_col_Idx;
	/*
	c++ 11:
	std::vector<cv::Mat> v_Mat = { m_maxPost_row_Idx, m_maxPost_col_Idx, m_down_row_Idx, m_down_col_Idx,
		m_recon_row_Idx, m_recon_col_Idx };*/
	int windowHeight = (r_idx_size + 10) * 2;// The height of the new composite image to be formed
	int nRows = 2; // Number of rows of images. (Number of columns will be calculated

	// depending on the value of total number of images).
	cv::Mat bigImg = makeCanvas(v_Mat, windowHeight, nRows);
	cv::namedWindow("Canvas");
	cv::imshow("Canvas", bigImg);
	cv::waitKey(0);
	}

			
	// thirdly, save the result
	// save image dimension
	v_dst_maxPost_row_column_Idx[imgIdx][0] = image_h; // image height
	v_dst_maxPost_row_column_Idx[imgIdx][1] = image_w; // image width
	// save indices
	for (int i = 0, tempIdx = 2; i < r_idx_size; ++i)
		for (int j = 0; j < c_idx_size; ++j, ++tempIdx)
			v_dst_maxPost_row_column_Idx[imgIdx][tempIdx] = m_recon_row_Idx.at<uchar>(i, j) * eWidth + m_recon_col_Idx.at<uchar>(i, j);

}


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
	){
	//******************************************************
	// first, read the input indicex into vector variable;
	//******************************************************
	// input max-posterior indices
	int numInputImage = 0; // number of images;
	vector<vector<int>> v_src_maxPost_row_column_Idx = read_row_col_idx_form_YML(whatkindImgs, s_src_Idx, numInputImage);

	//******************************************************
	// then, do uniform quantization to the learned indices;
	//******************************************************
	// function: vector (const vector& x); or vector (const vector& x, const allocator_type& alloc);
	// copy constructor (and copying with allocator)
	// Constructs a container with a copy of each of the elements in x, in the same order.
	vector<vector<int>> v_dst_maxPost_row_column_Idx(v_src_maxPost_row_column_Idx); // output de-quantized max-posterior indices
	vector<vector<int>> v_quantized_maxPost_row_column_Idx(v_src_maxPost_row_column_Idx); // quantized max-posterior indices

	for (int imgIdx = 0; imgIdx != numInputImage; ++imgIdx){/*for each image*/
		int r_idx_size = 0, c_idx_size = 0,
			image_h = v_src_maxPost_row_column_Idx[imgIdx][0], // image height
			image_w = v_src_maxPost_row_column_Idx[imgIdx][1]; // image width
		bool verbose = false;
#ifdef _DEBUG
		verbose = true;
#endif
		// here we do not need the return value of this function
		// here we do not need the return value of this function
		vector<int> v_rIdx, v_cIdx;
		getRowColIdxSizeOfOneImg(patchSideLength, patchSpacing, image_h, image_w, verbose, r_idx_size, c_idx_size, v_rIdx, v_cIdx);
		vector<int> v_maxPost_row_Idx(r_idx_size * c_idx_size, 0), v_maxPost_col_Idx(r_idx_size *c_idx_size, 0);

		// get row indices
		for (int i = 0, tempIdx = 2; i < r_idx_size; ++i)
			for (int j = 0; j < c_idx_size; ++j, ++tempIdx)
				v_maxPost_row_Idx[i * c_idx_size + j] = v_src_maxPost_row_column_Idx[imgIdx][tempIdx] / eWidth;
		// get column indices
		for (int i = 0, tempIdx = 2; i < r_idx_size; ++i)
			for (int j = 0; j < c_idx_size; ++j, ++tempIdx)
				v_maxPost_col_Idx[i * c_idx_size + j] = v_src_maxPost_row_column_Idx[imgIdx][tempIdx] % eWidth;


		vector <int> v_quantized_maxPost_row_Idx(r_idx_size *c_idx_size), // quantized row indices;
			v_quantized_maxPost_col_Idx(r_idx_size *c_idx_size); // quantized column indices;
		if (IsUniformQuantize){ // do uniform quantization;
			// do quantization and de-quantization
			quan_dequantize_row_col_Idx(r_idx_size, c_idx_size, eWidth, QuantizationLevel, v_maxPost_row_Idx, v_maxPost_col_Idx,
				v_quantized_maxPost_row_Idx, v_quantized_maxPost_col_Idx);
			// save the quantized indices into the quantized vector;
			for (int i = 0, tempIdx = 2; i < r_idx_size; ++i)
				for (int j = 0; j < c_idx_size; ++j, ++tempIdx)
					v_quantized_maxPost_row_column_Idx[imgIdx][tempIdx]
					= v_quantized_maxPost_row_Idx[i * c_idx_size + j] * eWidth + v_quantized_maxPost_col_Idx[i * c_idx_size + j];
			// save the de-quantized indices into the de-quantized or named as reconstruction result;
			for (int i = 0, tempIdx = 2; i < r_idx_size; ++i)
				for (int j = 0; j < c_idx_size; ++j, ++tempIdx)
					v_dst_maxPost_row_column_Idx[imgIdx][tempIdx]
					= v_maxPost_row_Idx[i * c_idx_size + j] * eWidth + v_maxPost_col_Idx[i * c_idx_size + j];
		} /* end of do uniform quantization*/
		else { // no uniform quantization;
			std::cout << "No uniform quantization to row/column indices!\n";
		}
		

		if (IsDownSampled){/*do Down- or Up-sampling*/
			// for v_quantized_maxPost_row_column_Idx
			down_up_sampling_2_row_col_Idx(r_idx_size, c_idx_size, imgIdx, image_h, image_w, eWidth, IsShowingImgFlag, v_quantized_maxPost_row_Idx,
				v_quantized_maxPost_col_Idx, f_width, f_height, v_quantized_maxPost_row_column_Idx);
			// for v_dst_maxPost_row_column_Idx;
			down_up_sampling_2_row_col_Idx(r_idx_size, c_idx_size, imgIdx, image_h, image_w, eWidth, IsShowingImgFlag, v_maxPost_row_Idx,
				v_maxPost_col_Idx, f_width,f_height, v_dst_maxPost_row_column_Idx);
		} /* end of do Down- or Up-sampling*/
		else{/*No Down- or Up-sampling*/
			std::cout << "No Down- or Up-sampling to row/column indices!\n";
		}/*End of No Down- or Up-sampling*/

		// release memory
		vector <int>().swap(v_quantized_maxPost_row_Idx);
		vector <int>().swap(v_quantized_maxPost_col_Idx);
		vector<int>().swap(v_maxPost_row_Idx);
		vector<int>().swap(v_maxPost_col_Idx);

	}/*end of each image*/

	//******************************************************************
	// third, save the uniform quantized and/or down-sampled/up-sampled 
	// indices to the output YML file;
	//******************************************************************
	write_row_col_idx_2_YML(whatkindImgs, s_dst_Idx,       v_dst_maxPost_row_column_Idx,       numInputImage);
	write_row_col_idx_2_YML(whatkindImgs, s_quantized_Idx, v_quantized_maxPost_row_column_Idx, numInputImage);
	// release memory
	vector<vector<int>>().swap(v_src_maxPost_row_column_Idx);
	vector<vector<int>>().swap(v_dst_maxPost_row_column_Idx);
	vector<vector<int>>().swap(v_quantized_maxPost_row_column_Idx);

}

	void getRowColIdxSizeOfOneImg(const int & patchSideLengh, // = 8
	const int & patchSpacing, // = 4, or 8
	const int & image_h, // image height
	const int & image_w, // image width
	const bool & verbose,
	int & r_idx_size,
	int & c_idx_size,
	vector<int> & v_rIdx, 
	vector<int> & v_cIdx
	){

	//////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////
	// generating possible patches extracted out of the current input image
	// last possible location of patches in the input image x;

	int c_xEndIdx = image_w - patchSideLengh;
	int r_xEndIdx = image_h - patchSideLengh;


	// possible indices of each xPatch in currently read image
	// std::vector<int> v_rIdx(0);
	// std::vector<int> v_cIdx(0);

	// clear element firstly
	v_rIdx.clear();
	v_cIdx.clear();
	if (verbose){
		std::cout << "\nc_xEndIdx = " << c_xEndIdx << ", r_xEndIdx = " << r_xEndIdx << "\n";
	}

	for (int r_Idx = 0; r_Idx < r_xEndIdx; r_Idx += patchSpacing){
		v_rIdx.push_back(r_Idx);
	}

	for (int c_Idx = 0; c_Idx < c_xEndIdx; c_Idx += patchSpacing){
		v_cIdx.push_back(c_Idx);
	}

	// include the last possible location of patches in the input image x;
	v_rIdx.push_back(r_xEndIdx);
	v_cIdx.push_back(c_xEndIdx);
	r_idx_size = v_rIdx.size();
	c_idx_size = v_cIdx.size();
	int PatchNum = r_idx_size* c_idx_size;

	if (verbose){
		std::cout << "\n The size of v_rIdx is: " << r_idx_size
			<< "\n The size of v_cIdx is: " << c_idx_size
			<< "\n The number of possible patches is: " << PatchNum
			<< endl;
	}

	//tuple<vector<int>, vector<int>> tuple_vector_out;
	//get<0>(tuple_vector_out) = v_rIdx;
	//get<1>(tuple_vector_out) = v_cIdx;
	//return tuple_vector_out;
}

// 1/4 down-sample of mat  
// recon = interpolated  up-sample of the 1/4 down-sample result,
// that is , the dimension of recon is the same as that of src.
cv::Mat quarter_down_sample_mat(cv::Mat & src, cv::Mat & recon){
	int h = src.rows, w = src.cols;
	int h_o = h / 2 + 1, w_o = w / 2 + 1;
	Mat out(h_o, w_o, CV_8UC1);
	int flag_row = (h % 2 == 0) ? 1 : 0, flag_col = (w % 2 == 0) ? 1 : 0;

	// first, do not include the last row and last column of the down_sample_resutl "out"
	for (int i = 0, i_size = h_o - 1; i < i_size; ++i) {
		for (int j = 0, j_size = w_o - 1; j < j_size; ++j){
			int temp = i * w_o + j;
			out.at<uchar>(i, j) = src.at<uchar>(2 * i, 2 * j);
			recon.at<uchar>(2 * i, 2 * j) = src.at<uchar>(2 * i, 2 * j);
		}
	}

	// then, for the last row, without containing the last element (i.e., right down corner element)
	int last_row_idx_out = h_o - 1;
	int last_row_idx_src = 2 * last_row_idx_out - flag_row;
	for (int j = 0, j_size = w_o - 1; j < j_size; ++j){
		out.at<uchar>(last_row_idx_out, j) = src.at<uchar>(last_row_idx_src, 2 * j);
		recon.at<uchar>(last_row_idx_src, 2 * j) = src.at<uchar>(last_row_idx_src, 2 * j);
	}

	// then, for the last column, without containing the last element (i.e., right down corner element)
	int last_col_idx_out = w_o - 1;
	int last_col_idx_src = 2 * last_col_idx_out - flag_col;
	for (int i = 0, i_size = h_o - 1; i < i_size; ++i){
		out.at<uchar>(i, last_col_idx_out) = src.at<uchar>(2 * i, last_col_idx_src);
		recon.at<uchar>(2 * i, last_col_idx_src) = src.at<uchar>(2 * i, last_col_idx_src);
	}

	// at last, the right down corner element
	out.at<uchar>(last_row_idx_out, last_col_idx_out) = src.at<uchar>(last_row_idx_src, last_col_idx_src);
	recon.at<uchar>(last_row_idx_src, last_col_idx_src) = src.at<uchar>(last_row_idx_src, last_col_idx_src);
	return out;
	// to be continued ...
	// to be continued ...

}

cv::Mat col_half_down_sample_mat(cv::Mat & src, cv::Mat & recon){
	int h = src.rows, w = src.cols;
	int h_o = h, w_o = w / 2 + 1;
	Mat out(h_o, w_o, CV_8UC1);
	int flag_col = (w % 2 == 0) ? 1 : 0;

	// first, do not include the last column of the down_sample_resutl "out"
	for (int i = 0, i_size = h_o; i < i_size; ++i) {
		for (int j = 0, j_size = w_o - 1; j < j_size; ++j){
			out.at<uchar>(i, j) = src.at<uchar>(i, 2 * j);
			recon.at<uchar>(i, 2 * j) = src.at<uchar>(i, 2 * j);
		}
	}

	// then, for the last column
	int last_col_idx_out = w_o - 1;
	int last_col_idx_src = 2 * last_col_idx_out - flag_col;
	for (int i = 0, i_size = h_o; i < i_size; ++i){
		out.at<uchar>(i, last_col_idx_out) = src.at<uchar>(i, last_col_idx_src);
		recon.at<uchar>(i, last_col_idx_src) = src.at<uchar>(i, last_col_idx_src);
	}

	// then interpolation of recon
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w - 2; j += 2){
			recon.at<uchar>(i, j + 1) = (src.at<uchar>(i, j) + src.at<uchar>(i, j + 2)) / 2;
		}
	}
	for (int i = 0; i < h; ++i) {
		for (int j = 0; j < w - 2; j += 2){
			recon.at<uchar>(i, j + 1) = (src.at<uchar>(i, j) + src.at<uchar>(i, j + 2)) / 2;
		}
	}
	return out;
}


// change row_column indices into j2k image for a less compressed file size
uintmax_t saveRowColIdx2_J2KImg(
	const int & patchSideLength,
	const int & patchSpacing,
	const int & eWidth,
	const string & whatkindImgs,
	const string & JPEG2K_EXE_Base_Path,
	const string & category_RowColIdx, // e.g., "Edwin-less_4_0_RowCol_Idx";
	const string & path_YML_Idx, // e.g., = ".../Edwin-less_4_RowCol_Idx.yml"
	const string & J2K_Result_BasePath, //
	const bool & is7zCompression, // do or not do compression from bmp files to 7z files;
	uintmax_t & size_7z_file // if do 7z compression, then find its file size; otherwise set the value of 0;
	){

	//*********************************************************************
	// read row and column indices from the above saved ".yml" file.
	//*********************************************************************

	int numInputImage;
	vector<vector<int>> v_maxPost_row_column_Idx = read_row_col_idx_form_YML(whatkindImgs, path_YML_Idx, numInputImage);
#ifdef _DEBUG
	cout << " numInputImage = " << numInputImage << endl;
#endif
	string JPEG2K_Compress_Path = JPEG2K_EXE_Base_Path + "/" + "opj_compress.exe";
	char * char_JPEG2K_Compress_Path = new char[JPEG2K_Compress_Path.length() + 1];
	std::strcpy(char_JPEG2K_Compress_Path, JPEG2K_Compress_Path.c_str());

	// make some necessary directories
	string tempDir = J2K_Result_BasePath + "/" + category_RowColIdx,
		tempDir_BMP = tempDir + "/" + "Idx_BMP",
		tempDir_J2K = tempDir + "/" + "Idx_J2K";
	MakeDir(tempDir);
	MakeDir(tempDir_BMP);
	MakeDir(tempDir_J2K);

	/*for each image*/
	for (int imgIdx = 0; imgIdx != numInputImage; ++imgIdx){
		int r_idx_size, c_idx_size,
			image_h = v_maxPost_row_column_Idx[imgIdx][0], // image height
			image_w = v_maxPost_row_column_Idx[imgIdx][1]; // image width
		bool verbose = false;
#ifdef _DEBUG
		verbose = true;
#endif
		// here we do not need the return value of this function
		vector<int> v_rIdx, v_cIdx;
		getRowColIdxSizeOfOneImg(patchSideLength, patchSpacing, image_h, image_w, verbose, r_idx_size, c_idx_size, v_rIdx, v_cIdx);
		vector<int> v_maxPost_row_Idx(r_idx_size * c_idx_size, 0), v_maxPost_col_Idx(r_idx_size *c_idx_size, 0);

		for (int i = 0, tempIdx = 2; i < r_idx_size; ++i)
			for (int j = 0; j < c_idx_size; ++j, ++tempIdx)
				v_maxPost_row_Idx[i * c_idx_size + j] = v_maxPost_row_column_Idx[imgIdx][tempIdx] / eWidth;


		for (int i = 0, tempIdx = 2; i < r_idx_size; ++i)
			for (int j = 0; j < c_idx_size; ++j, ++tempIdx)
				v_maxPost_col_Idx[i * c_idx_size + j] = v_maxPost_row_column_Idx[imgIdx][tempIdx] % eWidth;


		Mat m_maxPost_row_Idx(r_idx_size, c_idx_size, CV_8UC1), m_maxPost_col_Idx(r_idx_size, c_idx_size, CV_8UC1);


		for (int i = 0; i < r_idx_size; ++i)
			for (int j = 0; j < c_idx_size; ++j)
				m_maxPost_col_Idx.at<uchar>(i, j) = v_maxPost_col_Idx[i * c_idx_size + j];


		for (int i = 0; i < r_idx_size; ++i)
			for (int j = 0; j < c_idx_size; ++j)
				m_maxPost_row_Idx.at<uchar>(i, j) = v_maxPost_row_Idx[i * c_idx_size + j];
#ifdef _DEBUG
		double max_cIdx = .0, min_cIdx = .0,
			max_rIdx = .0, min_rIdx = .0;
		cv::minMaxLoc(m_maxPost_row_Idx, &min_rIdx, &max_rIdx);
		cv::minMaxLoc(m_maxPost_col_Idx, &min_cIdx, &max_cIdx);
		cout << "max_cIdx = " << max_cIdx
			<< " ,min_cIdx = " << min_cIdx
			<< ", max_rIdx = " << max_rIdx
			<< ", min_rIdx = " << min_rIdx << endl;
#endif


		char r_buffer[20], c_buffer[20];
		std::sprintf(r_buffer, "rowIdx_%03d", imgIdx);
		std::sprintf(c_buffer, "colIdx_%03d", imgIdx);


		// write idx to bmp images
		string rowbmpImg_path = tempDir_BMP + "/" + r_buffer + ".bmp",
			colbmpImg_path = tempDir_BMP + "/" + c_buffer + ".bmp",
			rowj2kImg_path = tempDir_J2K + "/" + r_buffer + ".j2k",
			colj2kImg_path = tempDir_J2K + "/" + c_buffer + ".j2k";
		cv::imwrite(rowbmpImg_path, m_maxPost_row_Idx);
		cv::imwrite(colbmpImg_path, m_maxPost_col_Idx);

		char * p_rowbmpImg_path = new char[rowbmpImg_path.length() + 1];
		char * p_colbmpImg_path = new char[colbmpImg_path.length() + 1];
		char * p_rowj2kImg_path = new char[rowj2kImg_path.length() + 1];
		char * p_colj2kImg_path = new char[colj2kImg_path.length() + 1];
		std::strcpy(p_rowbmpImg_path, rowbmpImg_path.c_str());
		std::strcpy(p_colbmpImg_path, colbmpImg_path.c_str());
		std::strcpy(p_rowj2kImg_path, rowj2kImg_path.c_str());
		std::strcpy(p_colj2kImg_path, colj2kImg_path.c_str());



		// losslessly encode bmp images into j2k images
		char * p_row_compress_buffer = new char[MAX_CHAR_NUM_OF_FILES_PATH];
		char * p_col_compress_buffer = new char[MAX_CHAR_NUM_OF_FILES_PATH];


		std::sprintf(p_row_compress_buffer, "%s -i %s -o %s", char_JPEG2K_Compress_Path, p_rowbmpImg_path, p_rowj2kImg_path);
		std::sprintf(p_col_compress_buffer, "%s -i %s -o %s", char_JPEG2K_Compress_Path, p_colbmpImg_path, p_colj2kImg_path);

		// compress, from .bmp to .j2k;
		std::system(p_row_compress_buffer);
		std::system(p_col_compress_buffer);

		delete[] p_rowbmpImg_path;
		delete[] p_colbmpImg_path;
		delete[] p_rowj2kImg_path;
		delete[] p_colj2kImg_path;
		vector <int>().swap(v_maxPost_row_Idx);
		vector <int>().swap(v_maxPost_col_Idx);

	} /*end of each-image-loop*/

	if (is7zCompression){
		// losslessly encode bmp images into .7z files
		string s_7z_Compress_Path = "C:/Users/ccj/Desktop/CPlusPlusProjects_CCJ/7z/7z.exe";
		char * char_7z_Compress_Path = new char[s_7z_Compress_Path.length() + 1];
		std::strcpy(char_7z_Compress_Path, s_7z_Compress_Path.c_str());
		char * input_bmp_dir = new char[tempDir_BMP.length() + 1];
		std::strcpy(input_bmp_dir, tempDir_BMP.c_str());
		string s_7z_output = tempDir_BMP + ".7z";
		char * output_7z_dir = new char[s_7z_output.length() + 1];
		std::strcpy(output_7z_dir, s_7z_output.c_str());

		char * p_7z_compress_buffer = new char[MAX_CHAR_NUM_OF_FILES_PATH];
		std::sprintf(p_7z_compress_buffer, "%s a %s %s", char_7z_Compress_Path, output_7z_dir, input_bmp_dir);

		// compress, from .bmp to .7z;
		std::system(p_7z_compress_buffer);
		// get file size
		size_7z_file = bf::file_size(s_7z_output);
		delete[] char_7z_Compress_Path;
		delete[] input_bmp_dir;
		delete[] output_7z_dir;
		delete[] p_7z_compress_buffer;

	}
	else {
		std::cout << "Do not do 7z compression to bmp files!\n";
		size_7z_file = 0;
	}

	vector<vector<int>>().swap(v_maxPost_row_column_Idx);
	delete[] char_JPEG2K_Compress_Path;
	// get j2k images file size
	return getFileSizeFromDir(tempDir_J2K);
}