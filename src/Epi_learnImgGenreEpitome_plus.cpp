#include "ImageEpitome.hpp"
#include "omp.h"

/*
Usage:
//#include <cstdlib>
char buffer[50];
int epitomeWidth = 256, patchSideLength = 8, patchSpacing = 4, numIteration = 10;
sprintf(buffer, "%03d_%02d_%02d_%02d-baseline", epitomeWidth, patchSideLength, patchSpacing, numIteration);
// buffer == "256_08_04_10-baseline"
*/

void ImgEpitome::learnImgGenreEpitome_plus(
	const vector<string> & v_imageFileNameList // read input images via their file-names.
	){

	// print the consuming time
	clock_t t = clock();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// some dimension parameters /////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// the total pixel number in the epitome
	int eNum = eHeight*eWidth;
	// the total pixel number in the image patch;
	int pNum = patchSideLengh * patchSideLengh;
	// the total number of the images in the current image category
	int ImageNum = v_imageFileNameList.size();

	char c_numIteration[20];
	sprintf(c_numIteration, "%02d", numIteration);
	std::string s_numIteration(c_numIteration);
	std::string fileStorageName = EpitomeResultDir + "/" + whatkindImgs + "_RowCol_Idx.yml";


	std::cout << "We are going to learn " << whatkindImgs << " category's Epitome."
		<< "\nThis category has " << ImageNum << " images in total.\n";

	// before the EM algorithm, keep the time
	t = clock() - t;
	std::printf("It took me %f seconds for initialization before EM.\n", ((float)t) / CLOCKS_PER_SEC);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// EM  Algorithm /////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// all the variables needed for the EM loop
	// Src(i.e., input signal) is Epitome, kernel signal is the patch extracted  out the original input images, or ones
	// *****************************
	// **** for Src and Kernel *****
	// *****************************


	c_fft::FFT_Workspace ws;
	c_fft::FFT_Convolution_Mode mode = c_fft::LINEAR_CORRELATION_VALIDSAME;
	c_fft::init_workspace(ws, mode, eHeight, eWidth, patchSideLengh, patchSideLengh, verbose);
	// since the following calculation might change the pointers ws.out_src and ws.out_kernel,
	// so when we call c_fft::clear_workspace(ws), the memorry address which will be deleted will not be found,
	// leading to error.
	// so we have to save the initial memory address, i.e., the value of pointers ws.out_src and ws.out_kernel.
	double * save_out_src;
	double * save_out_kernel;
	save_out_src = ws.out_src;
	save_out_kernel = ws.out_kernel;


	// shared variables initialization for the following parallel calculation;
	// Src Signal: epitome; kernel : image patches
	// constants
	v_patchSize_Ones = vector<Epitome_Doulbe>(pNum, 1.0); // ones
	if (isgrayScale){/*gray images*/
		/*Gray channel*/
		v_eMeanOverVar_Gray = vector<Epitome_Doulbe>(eNum);
		v_InvVar_Gray = vector<Epitome_Doulbe>(eNum);
		v_LogVar_Gray = vector<Epitome_Doulbe>(eNum); // not xPatch
		v_eMean2OverVar_Gray = vector<Epitome_Doulbe>(eNum); // not xPatch
		v_eLogVarSum_Gray = vector<Epitome_Doulbe>(eNum);
		v_eeSum_Gray = vector<Epitome_Doulbe>(eNum);
	}
	else {/*color images*/
		/*R channel*/
		v_eMeanOverVar_R = vector<Epitome_Doulbe>(eNum);
		v_InvVar_R = vector<Epitome_Doulbe>(eNum);
		v_LogVar_R = vector<Epitome_Doulbe>(eNum); // not xPatch
		v_eMean2OverVar_R = vector<Epitome_Doulbe>(eNum); // not xPatch
		v_eLogVarSum_R = vector<Epitome_Doulbe>(eNum);
		v_eeSum_R = vector<Epitome_Doulbe>(eNum);
		/*G channel*/
		v_eMeanOverVar_G = vector<Epitome_Doulbe>(eNum);
		v_InvVar_G = vector<Epitome_Doulbe>(eNum);
		v_LogVar_G = vector<Epitome_Doulbe>(eNum); // not xPatch
		v_eMean2OverVar_G = vector<Epitome_Doulbe>(eNum); // not xPatch
		v_eLogVarSum_G = vector<Epitome_Doulbe>(eNum);
		v_eeSum_G = vector<Epitome_Doulbe>(eNum);
		/*B channel*/
		v_eMeanOverVar_B = vector<Epitome_Doulbe>(eNum);
		v_InvVar_B = vector<Epitome_Doulbe>(eNum);
		v_LogVar_B = vector<Epitome_Doulbe>(eNum); // not xPatch
		v_eMean2OverVar_B = vector<Epitome_Doulbe>(eNum); // not xPatch
		v_eLogVarSum_B = vector<Epitome_Doulbe>(eNum);
		v_eeSum_B = vector<Epitome_Doulbe>(eNum);
	}

	// due to the limit of memory space, 
	// we separate the total PatchNum patches, into several chunks,
	// and each chunk has up to MAX_PARALLEL_SIZE patches which will further
	// involve the parallel computation.
	// // v_p_sumQ_sumQX_RGB_etc[7]  = { 
	//  p_sumQ, // summation of 3 channels;
	//  p_sumQX_R, p_sumQXX_R, 
	//  p_sumQX_G, p_sumQXX_G,
	//  p_sumQX_B, p_sumQXX_B,
	//}; 
	vector<vector<Epitome_Doulbe * const>> v_v_p_sumQ_sumQX_RGB_etc;
	if (isgrayScale){ /*gray */
		v_v_p_sumQ_sumQX_RGB_etc = vector<vector<Epitome_Doulbe * const>>(MAX_PARALLEL_SIZE, vector<Epitome_Doulbe * const>(3));
		for (int id = 0; id < MAX_PARALLEL_SIZE; ++id) {
			v_v_p_sumQ_sumQX_RGB_etc[id][0] = new Epitome_Doulbe[eNum]; // for  p_sumQ, #0
			v_v_p_sumQ_sumQX_RGB_etc[id][1] = new Epitome_Doulbe[eNum]; // for  p_sumQX_Gray; #1
			v_v_p_sumQ_sumQX_RGB_etc[id][2] = new Epitome_Doulbe[eNum]; // for  p_sumQXX_Gray; #2
		}
		// immediately do the zero initialization to each element, 
		// since later we will do the "+=" operation to each element,
		// so their initialized value must be ZERO.
		for (int ePixel_idx = 0; ePixel_idx < eNum; ++ePixel_idx)
			for (int id = 0; id < MAX_PARALLEL_SIZE; ++id){
				v_v_p_sumQ_sumQX_RGB_etc[id][0][ePixel_idx] = .0;
				v_v_p_sumQ_sumQX_RGB_etc[id][1][ePixel_idx] = .0;
				v_v_p_sumQ_sumQX_RGB_etc[id][2][ePixel_idx] = .0;
			}
		std::cout << "Allocating " << MAX_PARALLEL_SIZE << " blocks of zeros for parallel computation finishes!\n";

	}
	else {  /*color */
		v_v_p_sumQ_sumQX_RGB_etc = vector<vector<Epitome_Doulbe * const>>(MAX_PARALLEL_SIZE, vector<Epitome_Doulbe * const>(7));
		for (int id = 0; id < MAX_PARALLEL_SIZE; ++id) {
			v_v_p_sumQ_sumQX_RGB_etc[id][0] = new Epitome_Doulbe[eNum]; // for p_sumQ, #0
			v_v_p_sumQ_sumQX_RGB_etc[id][1] = new Epitome_Doulbe[eNum]; // for p_sumQX_R; #1
			v_v_p_sumQ_sumQX_RGB_etc[id][2] = new Epitome_Doulbe[eNum]; // for p_sumQXX_R; #2
			v_v_p_sumQ_sumQX_RGB_etc[id][3] = new Epitome_Doulbe[eNum]; // for p_sumQX_G; #3
			v_v_p_sumQ_sumQX_RGB_etc[id][4] = new Epitome_Doulbe[eNum]; // for p_sumQXX_G; #4
			v_v_p_sumQ_sumQX_RGB_etc[id][5] = new Epitome_Doulbe[eNum]; // for p_sumQX_B; #5
			v_v_p_sumQ_sumQX_RGB_etc[id][6] = new Epitome_Doulbe[eNum]; // for p_sumQXX_B; #6

		}
		// immediately do the zero initialization to each element, 
		// since later we will do the "+=" operation to each element,
		// so their initialized value must be ZERO.
		for (int ePixel_idx = 0; ePixel_idx < eNum; ++ePixel_idx)
			for (int id = 0; id < MAX_PARALLEL_SIZE; ++id){
				v_v_p_sumQ_sumQX_RGB_etc[id][0][ePixel_idx] = .0;
				v_v_p_sumQ_sumQX_RGB_etc[id][1][ePixel_idx] = .0;
				v_v_p_sumQ_sumQX_RGB_etc[id][2][ePixel_idx] = .0;
				v_v_p_sumQ_sumQX_RGB_etc[id][3][ePixel_idx] = .0;
				v_v_p_sumQ_sumQX_RGB_etc[id][4][ePixel_idx] = .0;
				v_v_p_sumQ_sumQX_RGB_etc[id][5][ePixel_idx] = .0;
				v_v_p_sumQ_sumQX_RGB_etc[id][6][ePixel_idx] = .0;
			}
		std::cout << "Allocating " << MAX_PARALLEL_SIZE << " blocks of zeros for parallel computation finishes!\n";
	} /*end of color */

	// FFT results
	int Num_FFT_Complex_Outputs = ws.h_fftw * (ws.w_fftw / 2 + 1);

	//#endif

	// EM starts here
	long long TotalNum = 0;
	for (int em_iterator = 0; em_iterator < numIteration; ++em_iterator){ // Beginning of EM Algorithm Iteration

		if (verbose){
			std::cout << "For category: " << whatkindImgs << " ," << em_iterator + 1 << "-th EM Iteration starts now!\n";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// The definition of all the vector data and FFTW plans /////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// initializing eMean and eVar, etc.
		if (isgrayScale){/*Gray channel*/
			// for v_LogVar
			for (int i = 0; i != eNum; ++i)
				v_LogVar_Gray[i] = log(v_eVar_Gray[i]);

			// for v_InvVar
			for (int i = 0; i != eNum; ++i)
				v_InvVar_Gray[i] = 1.0 / v_eVar_Gray[i];

			// for v_eMeanOverVar
			for (int i = 0; i != eNum; ++i)
				v_eMeanOverVar_Gray[i] = v_eMean_Gray[i] / v_eVar_Gray[i];

			// for v_eMean2OverVar
			for (int i = 0; i != eNum; ++i)
				v_eMean2OverVar_Gray[i] = v_eMean_Gray[i] * v_eMean_Gray[i] / v_eVar_Gray[i];
		}
		
		else{ /*Color images*/
		/*Red channel*/
		// for v_LogVar
		for (int i = 0; i != eNum; ++i)
			v_LogVar_R[i] = log(v_eVar_Red[i]);
		// for v_InvVar
		for (int i = 0; i != eNum; ++i)
			v_InvVar_R[i] = 1.0 / v_eVar_Red[i];
		// for v_eMeanOverVar
		for (int i = 0; i != eNum; ++i)
			v_eMeanOverVar_R[i] = v_eMean_Red[i] / v_eVar_Red[i];
		// for v_eMean2OverVar
		for (int i = 0; i != eNum; ++i)
			v_eMean2OverVar_R[i] = v_eMean_Red[i] * v_eMean_Red[i] / v_eVar_Red[i];

		/*Green channel*/
		// for v_LogVar
		for (int i = 0; i != eNum; ++i)
			v_LogVar_G[i] = log(v_eVar_Gre[i]);
		// for v_InvVar
		for (int i = 0; i != eNum; ++i)
			v_InvVar_G[i] = 1.0 / v_eVar_Gre[i];
		// for v_eMeanOverVar
		for (int i = 0; i != eNum; ++i)
			v_eMeanOverVar_G[i] = v_eMean_Gre[i] / v_eVar_Gre[i];
		// for v_eMean2OverVar
		for (int i = 0; i != eNum; ++i)
			v_eMean2OverVar_G[i] = v_eMean_Gre[i] * v_eMean_Gre[i] / v_eVar_Gre[i];
		
		/*Blue channel*/
		// for v_LogVar
		for (int i = 0; i != eNum; ++i)
			v_LogVar_B[i] = log(v_eVar_Blu[i]);
		// for v_InvVar
		for (int i = 0; i != eNum; ++i)
			v_InvVar_B[i] = 1.0 / v_eVar_Blu[i];
		// for v_eMeanOverVar
		for (int i = 0; i != eNum; ++i)
			v_eMeanOverVar_B[i] = v_eMean_Blu[i] / v_eVar_Blu[i];
		// for v_eMean2OverVar
		for (int i = 0; i != eNum; ++i)
			v_eMean2OverVar_B[i] = v_eMean_Blu[i] * v_eMean_Blu[i] / v_eVar_Blu[i];
		}

		///////////////////////////////////////////////////
		// to calculate v_eLogVarSum /////////////
		///////////////////////////////////////////////////
		ws.mode = c_fft::LINEAR_CORRELATION_VALIDSAME;
		update_workspace_few(ws);
		if (isgrayScale){
			c_fft::convolve(ws, &v_LogVar_Gray[0], &v_patchSize_Ones[0]);
			memcpy(&v_eLogVarSum_Gray[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
		}
		else{
			// R, G, B channels
			c_fft::convolve(ws, &v_LogVar_R[0], &v_patchSize_Ones[0]);
			memcpy(&v_eLogVarSum_R[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			c_fft::convolve(ws, &v_LogVar_G[0], &v_patchSize_Ones[0]);
			memcpy(&v_eLogVarSum_G[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			c_fft::convolve(ws, &v_LogVar_B[0], &v_patchSize_Ones[0]);
			memcpy(&v_eLogVarSum_B[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
		}
		
		///////////////////////////////////////////////////
		// to calculate v_eeSum /////////////
		///////////////////////////////////////////////////
		if (isgrayScale){
			c_fft::convolve(ws, &v_eMean2OverVar_Gray[0], &v_patchSize_Ones[0]);
			memcpy(&v_eeSum_Gray[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
		}
		else{
			// R, G, B channels
			c_fft::convolve(ws, &v_eMean2OverVar_R[0], &v_patchSize_Ones[0]);
			memcpy(&v_eeSum_R[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			c_fft::convolve(ws, &v_eMean2OverVar_G[0], &v_patchSize_Ones[0]);
			memcpy(&v_eeSum_G[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			c_fft::convolve(ws, &v_eMean2OverVar_B[0], &v_patchSize_Ones[0]);
			memcpy(&v_eeSum_B[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
		}



		// important parameters for calculating the Epitome model
		// before the starting of extracted training patches loop, 
		// please make sure zero-initialized doubles to 
		// v_double_sumQ(eNum), v_double_sumQX(eNum), and v_double_sumQXX(eNum).
		// i.e., clear the matrices used for collecting sufficient statistics.
		// zero-initialized doubles
		// those values are the final results,
		// since we do parallel calculation for all the patches extracted from one image;
		// e.g., for image 1, thus the temporary summation of sumQ, sumQX, or sumQXX is
		// just from this image; but the real sumQ, sumQX (R, G, B, or gray) and 
		// sumQXX (R, G, B, or gray), are gained through the summation of all the images;
		// // important parameters for calculating the Epitome model
		if (isgrayScale){
			for (int i = 0; i != eNum; ++i)
				v_double_sumQ[i] = .0;
			for (int i = 0; i != eNum; ++i)
				v_double_sumQX_Gray[i] = .0;
			for (int i = 0; i != eNum; ++i)
				v_double_sumQXX_Gray[i] = .0;
		}
		else {
		for (int i = 0; i != eNum; ++i)
			v_double_sumQ[i] = .0;
		for (int i = 0; i != eNum; ++i)
			v_double_sumQX_R[i] = .0;
		for (int i = 0; i != eNum; ++i)
			v_double_sumQXX_R[i] = .0;
		for (int i = 0; i != eNum; ++i)
			v_double_sumQX_G[i] = .0;
		for (int i = 0; i != eNum; ++i)
			v_double_sumQXX_G[i] = .0;
		for (int i = 0; i != eNum; ++i)
			v_double_sumQX_B[i] = .0;
		for (int i = 0; i != eNum; ++i)
			v_double_sumQXX_B[i] = .0;
		}

		// pay attention here:
		// do the same zero-initialization as above;
		// when each new EM iteration begins;
		if (isgrayScale){ /*gray */
			// zero initialization to each element, 
			// since later we will do the "+=" operation to each element,
			// so their initialized value must be ZERO.
			for (int ePixel_idx = 0; ePixel_idx < eNum; ++ePixel_idx)
				for (int id = 0; id < MAX_PARALLEL_SIZE; ++id){
					v_v_p_sumQ_sumQX_RGB_etc[id][0][ePixel_idx] = .0;
					v_v_p_sumQ_sumQX_RGB_etc[id][1][ePixel_idx] = .0;
					v_v_p_sumQ_sumQX_RGB_etc[id][2][ePixel_idx] = .0;
				}
		}
		else {  /*color */
			// do the zero initialization to each element, 
			// since later we will do the "+=" operation to each element,
			// so their initialized value must be ZERO.
			for (int ePixel_idx = 0; ePixel_idx < eNum; ++ePixel_idx)
				for (int id = 0; id < MAX_PARALLEL_SIZE; ++id){
					v_v_p_sumQ_sumQX_RGB_etc[id][0][ePixel_idx] = .0;
					v_v_p_sumQ_sumQX_RGB_etc[id][1][ePixel_idx] = .0;
					v_v_p_sumQ_sumQX_RGB_etc[id][2][ePixel_idx] = .0;
					v_v_p_sumQ_sumQX_RGB_etc[id][3][ePixel_idx] = .0;
					v_v_p_sumQ_sumQX_RGB_etc[id][4][ePixel_idx] = .0;
					v_v_p_sumQ_sumQX_RGB_etc[id][5][ePixel_idx] = .0;
					v_v_p_sumQ_sumQX_RGB_etc[id][6][ePixel_idx] = .0;
				}
		}

		// for parallel computation
		// vector of c_fft::FFT_Workspace ws
		vector<c_fft::FFT_Workspace> v_ws(MAX_PARALLEL_SIZE);
		vector<double *> v_save_out_src(MAX_PARALLEL_SIZE, NULL);
		vector<double *> v_save_out_kernel(MAX_PARALLEL_SIZE, NULL);
		// initialization of them
		for (int ws_idx = 0; ws_idx < MAX_PARALLEL_SIZE; ++ws_idx){
			c_fft::FFT_Convolution_Mode mode = c_fft::LINEAR_CORRELATION_VALIDSAME;
			c_fft::init_workspace(v_ws[ws_idx], mode, eHeight, eWidth, patchSideLengh, patchSideLengh, verbose);
			// since the following calculation might change the pointers ws.out_src and ws.out_kernel,
			// so when we call c_fft::clear_workspace(ws), the memorry address which will be deleted will not be found,
			// leading to error.
			// so we have to save the initial memory address, i.e., the value of pointers ws.out_src and ws.out_kernel.
			v_save_out_src[ws_idx] = v_ws[ws_idx].out_src;
			v_save_out_kernel[ws_idx] = v_ws[ws_idx].out_kernel;
		}

		
		///////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////
		// for each image of the current image genre,    //////////
		// to extract all the possible training patches  //////////
		///////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////
		for (int current_image_index = 0; current_image_index != ImageNum; ++current_image_index){
			Mat mat_currentInputImage =
				imread(DatabaseDir + "/" + v_imageFileNameList[current_image_index], Epi_Image_Read_Flag); // read a gray level image
			int channels = mat_currentInputImage.channels();
			if (verbose){
				std::cout << "\n EM Iteration :" << em_iterator + 1 << "/" << numIteration
					<< "; the current input image with " << channels << " channels :"
					<< current_image_index + 1 << "/" << ImageNum << " has been read!\n";
			}
			// size of input image
			int image_w = mat_currentInputImage.cols;
			int image_h = mat_currentInputImage.rows;

			vector<int> v_size(2);
			v_size[0] = image_h;
			v_size[1] = image_w;
			// show error if no image is read;
			if (mat_currentInputImage.empty())
				cerr << " > No input image is being read!" << endl;

			/////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////
			// generating possible patches extracted out of the current input image /////
			// last possible location of patches in the input image x;
			int c_xEndIdx = image_w - patchSideLengh,
				r_xEndIdx = image_h - patchSideLengh;

			// to get possible indices of each xPatch in currently read image
			tuple<vector<int>, vector<int>> t_v_row_col_Idx = getRowColIdxOfOneImg(image_h, image_w);
			vector<int> v_rIdx = get<0>(t_v_row_col_Idx);
			vector<int> v_cIdx = get<1>(t_v_row_col_Idx);
			int patchWiggle = (patchSpacing / 2) > 1 ? (patchSpacing / 2) : 1;
			//	if (patchWiggle == 0) patchWiggle = 1; // avoid the case of patchSpacing < 2.

			int r_size = v_rIdx.size(), c_size = v_cIdx.size();
			// number of all the possible patches extracted from the current image;
			int PatchNum = r_size* c_size;

			/////////////////////////////////////////////////////////////////////////////////////
			// for each training patch to generate several parameters        //////////////////
			////////////////////////////////////////////////////////////////////////////////////


			// generate the row index and column index for extracting the patch from the input image.
			// vector::end(), Returns an iterator referring to the past-the-end element in the vector container.
			// The past-the-end element is the theoretical element that would follow the last element in the vector. 
			// It does not point to any element, and thus shall be dereferenced.


			// if random wiggling
			// here we do not consider all those four cases, what we do is: firstly to generate a random row_column indices, 
			// secondly to let temp_r/temp_c = 0 or r_xEndIdx/c_xEndIdx in case of index exceeding.

			//	Condition 1 : rowIdx = True AND columnIdx = True, i.e., random shifting can occur to both row indices and column indices.
			//	Condition 2 : rowIdx = True AND columnIdx = False, i.e., random shifting can only occur to row indices, but not to column indices.
			//	Condition 3 : rowIdx = False AND columnIdx = True, i.e., random shifting can only occur to column indices, but not to row indices.
			//	Condition 4 : rowIdx = False AND columnIdx = False, i.e., random shifting can not occur to neither row indices nor column indices.

			for (int r_it = 1; r_it < r_size - 1; r_it++) {
				// Here r_it != 0 AND r_it != r_size - 1, means that we randomize all the row indices except the first one and the last one. 
				// Specifically, if we randomize the regular indices starting with 0 (e.g., [0 4 8 ... 224 228 232 ]), we might not get zero-based sequences any more
				// (e.g., [ 2 5 9 ...222 229 232]), and this non-zero-based indices will affect the following image reconstruction.
				// Similarly for the last index of the row-indices.
				if (RandomShifting){
					int temp_r = v_rIdx[r_it] + rand() % (2 * patchWiggle + 1) - patchWiggle; // generate a random variable within [-patchWiggle, patchWiggle] 
					// the following trick is needed, in case of indices exceeding, 
					// especially when the last two elements of v_rIdx (or v_cIdx) are close, for example, they are 114 and 145, 
					// so if ( 144 + a random value within [-patchWiggle, patchWiggle]), 144 + 2 = 146 > 145, which will result in a case of indices exceeding.
					temp_r = temp_r < r_xEndIdx ? temp_r : r_xEndIdx;
					temp_r = temp_r >= 0 ? temp_r : 0;
					//save the random indices instead of the regular indices in v_rIdx and v_cIdx
					v_rIdx[r_it] = temp_r;
				}
			}

			for (int c_it = 1; c_it < c_size - 1; c_it++) { // Similarly for c_it != 0 AND c_it != c_size - 1.
				if (RandomShifting){
					int temp_c = v_cIdx[c_it] + rand() % (2 * patchWiggle + 1) - patchWiggle; // generate a random variable within [-patchWiggle, patchWiggle]
					// the following trick is needed, in case of indices exceeding, 
					// especially when the last two elements of v_rIdx (or v_cIdx) are close, for example, they are 114 and 145, 
					// so if ( 144 + a random value within [-patchWiggle, patchWiggle]), 144 + 2 = 146 > 145, which will result in a case of indices exceeding.
					temp_c = temp_c < c_xEndIdx ? temp_c : c_xEndIdx;
					temp_c = temp_c >= 0 ? temp_c : 0;
					//save the random indices instead of the regular indices in v_rIdx and v_cIdx
					v_cIdx[c_it] = temp_c;
				}
			}

			// the equivalent 1-D indices calculated from the 2-D coordiantes;
			// that is, temp_location_Idx = (temp_r, temp_c) =  temp_r * image_w + temp_c;
			// so, temp_r = temp_location_Idx / image_w ;
			// temp_c = temp_location_Idx % image_w;
			vector <int> v_1D_r_c_Idx(PatchNum);
			for (int r_it = 0, temp = 0; r_it < r_size; ++r_it) {
				for (int c_it = 0; c_it < c_size; ++c_it, ++ temp){
					//v_1D_r_c_Idx[r_it * c_size + c_it] = v_rIdx[r_it] * image_w + v_cIdx[c_it];
					v_1D_r_c_Idx[temp] = v_rIdx[r_it] * image_w + v_cIdx[c_it];
				}
			}

			int max_1D_r_c_Idx = r_xEndIdx *image_w + c_xEndIdx;


			/*
			// save the vector eMean into txt file
			if (em_iterator == 0){
				fstream fs_1D_Idx;
				string  binFileName = EpitomeResultDir + "/" + whatkindImgs + "_FFT_1D_Idx_" + s_numIteration + ".txt";
				fs_1D_Idx.open(binFileName, std::fstream::binary | std::fstream::out);

				if (fs_1D_Idx.is_open()){
					for (int r_it = 0; r_it < 10; r_it++){
						for (int c_it = 0; c_it < 10; c_it++){
							int i = r_it * c_size + c_it;
							fs_1D_Idx << v_1D_r_c_Idx[i] << " , " << " r = " << v_1D_r_c_Idx[i] / image_w << " ,  c = " << v_1D_r_c_Idx[i] % image_w  << ", ";
						}
						fs_1D_Idx << std::endl;
					}

					// Note: For binary files, reading and writing data with the extraction and insertion operators (<< and >>) and functions like getline is not efficient, 
					// since we do not need to format any data and data is likely not formatted in lines.
					//	 fs<< v_maxPost_row_column_Idx[i][j] << " ";
				}
				else cout << "Unable to open file";
				fs_1D_Idx.close();
			}
			*/

			// due to the limit of memory space, 
			// we separate the total PatchNum patches, into several chunks,
			// and each chunk has up to MAX_PARALLEL_SIZE patches which will further
			// involve the parallel computation.
			int total_paral_times = std::ceil((double)PatchNum / MAX_PARALLEL_SIZE);
			// but we have to add extra elements to "patchNum", so that it can be divided by MAX_PARALLEL_SIZE
			// with no remainder.
			PatchNum = total_paral_times * MAX_PARALLEL_SIZE;
			while (PatchNum > v_1D_r_c_Idx.size()){
				// rand() % max_1D_r_c_Idx, generate a random number, which is :
				// 0<< rand() % max_1D_r_c_Idx < max_1D_r_c_Idx;
				// v_1D_r_c_Idx.push_back((rand() % max_1D_r_c_Idx));
				// this method only guarantee the 1D_r_c_Idx <= max_1D_r_c_Idx,
				// but, after the 1-D indices changed to 2-D row/column,
				// that is, row = 1D_r_c_Idx / image_w , column =  1D_r_c_Idx % image_w;
				// the calculated row and column will exceed their limited value -- c_xEndIdx and  r_xEndIdx;
				// so we have to directly find the 2-D indices what we need, and then change them to 1-d indices;
				/*
				int rand_c = rand() % (c_xEndIdx + 1);
				int rand_r = rand() % (r_xEndIdx + 1);
				int temp_idx = rand_r * image_w + rand_c;
				v_1D_r_c_Idx.push_back(temp_idx);
				*/
				v_1D_r_c_Idx.push_back((rand() % (c_xEndIdx + 1)) + (rand() % (r_xEndIdx + 1)) * image_w);
			}
#ifdef _DEBUG
			std::cout << "before random value added, max_1D_r_c_Idx = " << max_1D_r_c_Idx << endl;
			std::cout << "After random value added, max_1D_r_c_Idx = " << *max_element(v_1D_r_c_Idx.begin(), v_1D_r_c_Idx.end()) << endl;
#endif // DEBUG

			int nThreads = 1;

			// the corresponding row_column indices of the maximum posterior element 
			// for all the possible xPatch extracted out of the current image
			if (em_iterator == (numIteration - 1)){
				// for the last EM iteration
				// we will save the maxPost_row_column_Idx
				v_maxPost_row_column_Idx[current_image_index] = vector<int>(PatchNum + 2, 0);
				v_maxPost_row_column_Idx[current_image_index][0] = image_h;
				v_maxPost_row_column_Idx[current_image_index][1] = image_w;
			}

			// for saving the max_unnormalized_post for each image patch out of one image;
			v_max_unnormalized_post = vector<Epitome_Doulbe>(PatchNum);

			//****************************************************************************
			//****************************************************************************
			//parallel calculation begins here for each xPatch(R, G, B, and/or Gray)
			//****************************************************************************
			//****************************************************************************

			long int trainingCounter = 0; // counting the number of training patches for showing the remaining time;
			

			// const bool IsSaveRowColIdx = (em_iterator == (numIteration - 1)) ? true : false;
			bool IsSaveRowColIdx = false;
			// due to the limit of memory space, 
			// we separate the total PatchNum patches, into several chunks,
			// and each chunk has up to MAX_PARALLEL_SIZE patches which will further
			// involve the parallel computation.

			
			omp_set_num_threads(NUM_THREADS);
			/*several parallel computation*/
			for (int paral_time = 0; paral_time < total_paral_times; ++paral_time){ 
				// par attention to those calculation,
				// _1D_Idx_iter = idx_start; _1D_Idx_iter < idx_end; ++_1D_Idx_iter
				// they can guarantee that _1D_Idx_iter can be every number in the range [0, PatchNum-1]
				int idx_start = paral_time* MAX_PARALLEL_SIZE;
				int idx_end = (paral_time + 1)* MAX_PARALLEL_SIZE;
				
				int _1D_Idx_iter;

			// omp_set_num_threads(NUM_THREADS);
			//omp_set_num_threads(1);

#pragma omp parallel  default(none) private (_1D_Idx_iter) shared (v_ws, idx_start, idx_end, nThreads, trainingCounter, v_1D_r_c_Idx, IsSaveRowColIdx, current_image_index, PatchNum, ImageNum, image_w, mat_currentInputImage,v_v_p_sumQ_sumQX_RGB_etc)
			{


#ifdef _DEBUG
#pragma omp master
				nThreads = omp_get_num_threads();
				if (nThreads == NUM_THREADS) {
#pragma omp master
					printf_s("%d OpenMP threads were used.\n", NUM_THREADS);
				}
				else {
#pragma omp master
					printf_s("Expected %d OpenMP threads, but %d were used.\n",
						NUM_THREADS, nThreads);
				}
#endif

#pragma omp for 
				for (_1D_Idx_iter = idx_start; _1D_Idx_iter < idx_end; ++_1D_Idx_iter){
					// paral_time
					
					bool flag = (_1D_Idx_iter == 998) ? true : false;
					int temp = 0;
					if (get_xPatch_sumQXX_etc_RGB_or_Gray_Channel(
						v_ws[_1D_Idx_iter - idx_start],
						_1D_Idx_iter, // used saving maxPost_row_column_Idx, for identifying each patch idx;
						// to save the maxPost_row_column_Idx or not, 
						// since we just want to save the indices during the last time EM iteration;
						IsSaveRowColIdx,
						v_1D_r_c_Idx[_1D_Idx_iter], //  temp_location = temp_r * image_w + temp_c;
						image_w,
						current_image_index, // used for saving maxPost_row_column_Idx;
						mat_currentInputImage,
						PatchNum, // used for disply the remaining time; 
						&v_max_unnormalized_post[_1D_Idx_iter],
						// v_p_sumQ_sumQX_RGB_etc[7]  = { 
						//  p_sumQ, // summation of 3 channels;
						//  p_sumQX_R, p_sumQXX_R, 
						//  p_sumQX_G, p_sumQXX_G,
						//  p_sumQX_B, p_sumQXX_B,
						//}; 

						// since v_v_p_sumQ_sumQX_RGB_etc have only up to MAX_PARALLEL_SIZE elements;
						// so pay attention to its indices;
						v_v_p_sumQ_sumQX_RGB_etc[_1D_Idx_iter - idx_start]
						))
#pragma omp critical 
						trainingCounter++; // if true, execute trainingCounter++
					////////////////////////////////////////////////////////////////////////////////////
					//************************************************************************************
					//	display and estimate the remaining time to finish the process of learning an Epitome model
					//************************************************************************************
					////////////////////////////////////////////////////////////////////////////////////
					// 
					
				}/*end of omp for-loop, each xPatch*/
#ifdef _DEBUG
				if ((verbose) && (trainingCounter % (PatchNum / 10) == 0))	
#pragma omp critical
					printf_s("For image %d of %d 's patch learning, %f %% Complete...\n", current_image_index, 
					ImageNum, 100 * trainingCounter / PatchNum);
#endif
			} /*end of omp parallel */

			

			// change to serial computation, i.e., 1 thread;
			//omp_set_num_threads(1);
			// and check whether it has been successfully set;
			//nThreads = omp_get_num_threads();
#ifdef _DEBUG
#pragma omp master
			printf_s("Parallel calculation finishes. Now %d OpenMP threads were used.\n", nThreads);
#endif
			
			}/*end of several parallel computation*/
			

			nThreads = omp_get_num_threads();
			// change to serial computation, i.e., 1 thread;
			if (nThreads != 1)
				omp_set_num_threads(1);

			// release some variable for next image operation
			// delete possible indices of each xPatch in currently read image
			std::vector<int>().swap(v_rIdx);
			std::vector<int>().swap(v_cIdx);

			TotalNum += PatchNum;

		} // end of each current input image



		// release vector of c_fft::FFT_Workspace ws
		for (int ws_idx = 0; ws_idx < MAX_PARALLEL_SIZE; ++ws_idx){
			// clear ws
			// get the initial pointer value
			v_ws[ws_idx].out_src = v_save_out_src[ws_idx];
			v_ws[ws_idx].out_kernel = v_save_out_kernel[ws_idx];
			c_fft::clear_workspace(v_ws[ws_idx]);
		}
		vector<c_fft::FFT_Workspace>().swap(v_ws);

		// accumulation for sum, sumQX, sumQXX, etc
		// v_v_p_sumQ_sumQX_RGB_etc
		// accumulate 

		for (int i = 0; i != eNum; ++i){
			for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
				v_double_sumQ[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][0][i];
			}
		}
		if (isgrayScale){
			for (int i = 0; i != eNum; ++i){
				for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
					v_double_sumQX_Gray[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][1][i];
				}
			}

			for (int i = 0; i != eNum; ++i){
				for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
					v_double_sumQXX_Gray[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][2][i];
				}
			}
		}
		else { /*color images*/
			for (int i = 0; i != eNum; ++i){
				for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
					v_double_sumQX_R[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][1][i];
				}
			}

			for (int i = 0; i != eNum; ++i){
				for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
					v_double_sumQXX_R[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][2][i];
				}
			}
			for (int i = 0; i != eNum; ++i){
				for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
					v_double_sumQX_G[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][3][i];
				}
			}

			for (int i = 0; i != eNum; ++i){
				for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
					v_double_sumQXX_G[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][4][i];
				}
			}

			for (int i = 0; i != eNum; ++i){
				for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
					v_double_sumQX_B[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][5][i];
				}
			}

			for (int i = 0; i != eNum; ++i){
				for (int paralPatchIdx = 0; paralPatchIdx < MAX_PARALLEL_SIZE; ++paralPatchIdx){
					v_double_sumQXX_B[i] += v_v_p_sumQ_sumQX_RGB_etc[paralPatchIdx][6][i];
				}
			}

		}/*end of color images*/

		// avoid numerical problems
		// to find the nearly -being-zero elements in sumQ
		for (int i = 0; i != eNum; ++i){
			if (v_double_sumQ[i] <= tolerance){
				v_double_sumQ[i] = 1;
				if (isgrayScale){
					v_double_sumQX_Gray[i] = v_eMean_Gray[i];
					v_double_sumQXX_Gray[i] = v_eVar_Gray[i] + v_eMean_Gray[i] * v_eMean_Gray[i];
				}
				else{
				v_double_sumQX_R[i] = v_eMean_Red[i];
				v_double_sumQXX_R[i] = v_eVar_Red[i] + v_eMean_Red[i] * v_eMean_Red[i];
				v_double_sumQX_G[i] = v_eMean_Gre[i];
				v_double_sumQXX_G[i] = v_eVar_Gre[i] + v_eMean_Gre[i] * v_eMean_Gre[i];
				v_double_sumQX_B[i] = v_eMean_Blu[i];
				v_double_sumQXX_B[i] = v_eVar_Blu[i] + v_eMean_Blu[i] * v_eMean_Blu[i];
				}
			}
		}

#ifdef _DEBUG
		int temp = em_iterator;
		std::cout << "Current EM " << em_iterator + 1 << "/ " << numIteration << endl;
#endif // _DEBUG

		// compute the new epitome
		if (isgrayScale){/*gray images*/
			for (int i = 0; i != eNum; ++i){
				v_eMean_Gray[i] = v_double_sumQX_Gray[i] / v_double_sumQ[i];
				v_eVar_Gray[i] = v_double_sumQXX_Gray[i] / v_double_sumQ[i] - v_eMean_Gray[i] * v_eMean_Gray[i];
				
				// eMean belongs to [0, 1] interval.
				// make sure that the mean is within the range 0 - 1 (may not because of numerical issues)
				//	eMean[i] = eMean[i] > 0 ? eMean[i] : 0;
				//	eMean[i] = eMean[i] < 1 ? eMean[i] : 1;
				v_eMean_Gray[i] = v_eMean_Gray[i] > 0 ? v_eMean_Gray[i] : 0;
				v_eMean_Gray[i] = v_eMean_Gray[i] < 1 ? v_eMean_Gray[i] : 1;

				//enforce a minimum variance
				//	eVar[i] = eVar[i] > minVar ? eVar[i] : minVar;
				v_eVar_Gray[i] = v_eVar_Gray[i] > minVar ? v_eVar_Gray[i] : minVar;
			}
		}/*end of gray images*/
		else {/*color images*/
		for (int i = 0; i != eNum; ++i){
			v_eMean_Red[i] = v_double_sumQX_R[i] / v_double_sumQ[i];
			v_eVar_Red[i] = v_double_sumQXX_R[i] / v_double_sumQ[i] - v_eMean_Red[i] * v_eMean_Red[i];
			v_eMean_Gre[i] = v_double_sumQX_G[i] / v_double_sumQ[i];
			v_eVar_Gre[i] = v_double_sumQXX_G[i] / v_double_sumQ[i] - v_eMean_Gre[i] * v_eMean_Gre[i];
			v_eMean_Blu[i] = v_double_sumQX_B[i] / v_double_sumQ[i];
			v_eVar_Blu[i] = v_double_sumQXX_B[i] / v_double_sumQ[i] - v_eMean_Blu[i] * v_eMean_Blu[i];

			// eMean belongs to [0, 1] interval.
			// make sure that the mean is within the range 0 - 1 (may not because of numerical issues)
			//	eMean[i] = eMean[i] > 0 ? eMean[i] : 0;
			//	eMean[i] = eMean[i] < 1 ? eMean[i] : 1;
			v_eMean_Red[i] = v_eMean_Red[i] > 0 ? v_eMean_Red[i] : 0;
			v_eMean_Red[i] = v_eMean_Red[i] < 1 ? v_eMean_Red[i] : 1;
			v_eMean_Gre[i] = v_eMean_Gre[i] > 0 ? v_eMean_Gre[i] : 0;
			v_eMean_Gre[i] = v_eMean_Gre[i] < 1 ? v_eMean_Gre[i] : 1;
			v_eMean_Blu[i] = v_eMean_Blu[i] > 0 ? v_eMean_Blu[i] : 0;
			v_eMean_Blu[i] = v_eMean_Blu[i] < 1 ? v_eMean_Blu[i] : 1;

			//enforce a minimum variance
			//	eVar[i] = eVar[i] > minVar ? eVar[i] : minVar;
			v_eVar_Red[i] = v_eVar_Red[i] > minVar ? v_eVar_Red[i] : minVar;
			v_eVar_Gre[i] = v_eVar_Gre[i] > minVar ? v_eVar_Gre[i] : minVar;
			v_eVar_Blu[i] = v_eVar_Blu[i] > minVar ? v_eVar_Blu[i] : minVar;
		}
		}/*end of color images*/


		std::cout << "\nThe " << em_iterator + 1 << "-th EM Iteration finishes now!\n";

	} // End of EM Algorithm

#ifdef _DEBUG
	std::cout << " TotalNum = " << TotalNum << endl;
#endif
	

	// after the EM algorithm, keep the time
	t = clock() - t;
	double temp_t = ((float)t) / CLOCKS_PER_SEC;
	printf_s("EM Iteration took %f seconds.\n", temp_t);
	
	fstream fs_time;
	string  timeFileName = EpitomeResultDir + "/" + whatkindImgs + "_timeRecord" + s_numIteration + ".txt";
	fs_time.open(timeFileName, std::fstream::out | std::fstream::app);
	if (fs_time.is_open()){
		fs_time
			<< " learning epitome in this experiment takes " << temp_t << " seconds." << endl
			<< " This experiment's parameters include: " << endl
			<< "   MAX_PARALLEL_SIZE = " << MAX_PARALLEL_SIZE << endl
			<< "   NUM_THREADS = " << NUM_THREADS << endl
			<< "   numIteration = " << numIteration << endl
			<< "   ImageNum = " << ImageNum << endl
			<< "   eHeight = " << eHeight << endl;
		if (isgrayScale)
		fs_time
			<< "   IsgrayScale = true, i.e., we have learned gray-scale images." << endl;
		else 
		fs_time
		    << "   IsgrayScale = false, i.e., we have learned color images." << endl;
		
		fs_time << "   TotalNum of image patches = " << TotalNum << endl << endl << endl 
			;
	}

	fs_time.close();

	// release containers of private variables for parallel computation
	for (int id = 0; id < MAX_PARALLEL_SIZE; ++id) {
		int j_size = v_v_p_sumQ_sumQX_RGB_etc[id].size();
		for (int j = 0; j < j_size; ++j)
			delete[] v_v_p_sumQ_sumQX_RGB_etc[id][j];
	}

	// clear ws
	// get the initial value
	ws.out_src = save_out_src;
	ws.out_kernel = save_out_kernel;
	c_fft::clear_workspace(ws);

	// memory release time
	t = clock() - t;
	printf_s("It took %f seconds for memory release.\n", ((float)t) / CLOCKS_PER_SEC);


	if (isgrayScale){ /*gray epitome*/
	Mat_<double> mat_eMean(eHeight, eWidth);
	
	// save the vector eMean and eVar into to mat_eMean and mat_eVar for the following txtFile saving
	for (int r = 0; r < eHeight; r++){
		for (int c = 0; c < eWidth; c++){
			// using operator [] to access element of vector
			int temp = r*eWidth + c;
			mat_eMean.at<double>(r, c) = v_eMean_Gray[temp];
		}
	}

	// save mat_eMean into to an image
	cv::imwrite(EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMean_EM" + s_numIteration + ".bmp", 255 * mat_eMean);
	
	vector<int> compression_params_png(2);
	compression_params_png[0] = CV_IMWRITE_PNG_COMPRESSION;
	compression_params_png[1] =0;
	cv::imwrite(EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMean_EM" + s_numIteration + ".png", 255 * mat_eMean, compression_params_png);

	vector<int> compression_params_jpeg(2);
	compression_params_jpeg[0] = CV_IMWRITE_JPEG_QUALITY;
	compression_params_jpeg[1] = 100;
	cv::imwrite(EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMean_EM" + s_numIteration + ".jpg", 255 * mat_eMean, compression_params_jpeg);

	// save the vector eMean into binary file
	fstream fs_eMean;
	string  binFileName = EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMean_EM" + s_numIteration + ".bin";
	fs_eMean.open(binFileName, std::fstream::binary | std::fstream::out);

	if (fs_eMean.is_open()){
		fs_eMean.write(reinterpret_cast<char*>(&v_eMean_Gray[0]), eHeight*eWidth*sizeof(Epitome_Doulbe));
		// Note: For binary files, reading and writing data with the extraction and insertion operators (<< and >>) 
		// and functions like getline is not efficient, 
		// since we do not need to format any data and data is likely not formatted in lines.
		//	 fs<< v_maxPost_row_column_Idx[i][j] << " ";
	}
	else std::cout << "Unable to open file";
	fs_eMean.close();

	// save the vector eMean and eVar into ".yml" file for image reconstruction
	FileStorage fs_Epitome_in(EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMeanVar_EM" + s_numIteration + ".yml", FileStorage::WRITE);
	if (!fs_Epitome_in.isOpened()){
		cerr << "failed to open " << EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMeanVar_EM" + s_numIteration + ".yml" << endl;
		fileStorageHelp();
	}
	else
		fs_Epitome_in << "e_Mean_Gray" << v_eMean_Gray << "e_Var_Gray" << v_eVar_Gray;
	fs_Epitome_in.release();

	if (verbose){
		std::cout
			<< " The newly learned v_eMean and v_eVar has been saved into "
			+ whatkindImgs + "_Gray_eMeanVar_EM" + s_numIteration + ".yml\n"
			<< " The newly learned v_eMean and v_eVar has been saved into "
			+ whatkindImgs + "_Gray_eMean_EM" + s_numIteration + ".bin\n"
			<< " The newly learned v_eMean has been saved as a JPG, BMP, PNG image.\n"
			<< " EM has finished! TotalNum = " << TotalNum << " patches.\n";
	}

	}/*end of gray epitome*/
	else { /*color epitome*/

		Mat mat_eMean(eHeight, eWidth, CV_64FC3);
		Mat mat_eVar(eHeight, eWidth, CV_64FC3);

		// save the vector eMean and eVar into to mat_eMean and mat_eVar for the following txtFile saving
		for (int r = 0; r < eHeight; r++){
			for (int c = 0; c < eWidth; c++){
				// using operator [] to access element of vector
				int temp = r*eWidth + c;
				mat_eMean.at<cv::Vec3d>(r, c)[0] = v_eMean_Blu[temp]; // blue
				mat_eMean.at<cv::Vec3d>(r, c)[1] = v_eMean_Gre[temp]; // green
				mat_eMean.at<cv::Vec3d>(r, c)[2] = v_eMean_Red[temp]; // red
				mat_eVar.at<cv::Vec3d>(r, c)[0] = v_eVar_Blu[temp]; // blue
				mat_eVar.at<cv::Vec3d>(r, c)[1] = v_eVar_Gre[temp]; // green
				mat_eVar.at<cv::Vec3d>(r, c)[2] = v_eVar_Red[temp]; //  red
			}
		}

		// save mat_eMean into to an image
		cv::imwrite(EpitomeResultDir + "/" + whatkindImgs + "_Color_eMean_EM" + s_numIteration + ".bmp", 255 * mat_eMean);
		vector<int> compression_params_png(2);
	    compression_params_png[0] = CV_IMWRITE_PNG_COMPRESSION;
	    compression_params_png[1] =0;
		cv::imwrite(EpitomeResultDir + "/" + whatkindImgs + "_Color_eMean_EM" + s_numIteration + ".png", 255 * mat_eMean, compression_params_png);

		vector<int> compression_params_jpeg(2);
	    compression_params_jpeg[0] = CV_IMWRITE_JPEG_QUALITY;
	    compression_params_jpeg[1] = 100;
		cv::imwrite(EpitomeResultDir + "/" + whatkindImgs + "_Color_eMean_EM" + s_numIteration + ".jpg", 255 * mat_eMean, compression_params_jpeg);

		// save the vector eMean into binary file
		fstream fs_eMean;
		string  binFileName = EpitomeResultDir + "/" + whatkindImgs + "_Color_eMean_EM" + s_numIteration + ".bin";
		fs_eMean.open(binFileName, std::fstream::binary | std::fstream::out);

		if (fs_eMean.is_open()){
			fs_eMean.write(reinterpret_cast<char*>(mat_eMean.data), mat_eMean.channels() * eHeight*eWidth*sizeof(Epitome_Doulbe));
			// Note: For binary files, reading and writing data with the extraction and insertion operators (<< and >>) and functions like getline is not efficient, 
			// since we do not need to format any data and data is likely not formatted in lines.
			//	 fs<< v_maxPost_row_column_Idx[i][j] << " ";
		}
		else std::cout << "Unable to open file";
		fs_eMean.close();

		// save the vector eMean and eVar into ".yml" file for image reconstruction
		FileStorage fs_Epitome_in(EpitomeResultDir + "/" + whatkindImgs + "_Color_eMeanVar_EM" + s_numIteration + ".yml", FileStorage::WRITE);
		if (!fs_Epitome_in.isOpened()){
			cerr << "failed to open " << EpitomeResultDir + "/" + whatkindImgs + "_Color_eMeanVar_EM" + s_numIteration + ".yml" << endl;
			fileStorageHelp();
		}
		else
			fs_Epitome_in << "e_Mean_Red" << v_eMean_Red << "e_Var_Red" << v_eVar_Red
			<< "e_Mean_Gre" << v_eMean_Gre << "e_Var_Gre" << v_eVar_Gre
			<< "e_Mean_Blu" << v_eMean_Blu << "e_Var_Blu" << v_eVar_Blu;
		fs_Epitome_in.release();

		if (verbose){
			std::cout
				<< " The newly learned v_eMean and v_eVar has been saved into " 
				+ whatkindImgs + "_Color_eMeanVar_EM" + s_numIteration + ".yml\n"
				<< " The newly learned v_eMean and v_eVar has been saved into "
				+ whatkindImgs + "_Color_eMean_EM" + s_numIteration + ".bin\n"
				<< " The newly learned v_eMean has been saved as a JPG, BMP, PNG image.\n"
				<< " EM has finished! TotalNum = " << TotalNum << " patches.\n";
		}
	} /*end of saving color epitome*/

	// to be continued .......
	// to be continued .......
	// to be continued .......

}