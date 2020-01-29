#include "ImageEpitome.hpp"

/*
std::string DatabaseDir,
std::string EpitomeResultDir,
std::string ReconsCompresDir,
std::string whatkindImgs,
*/
void ImgEpitome::learnImgGenreEpitome(
	const vector<string> & v_imageFileNameList // read input images via their file-names.
	){

	// print the consuming time
	clock_t t = clock();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// some dimension parameters /////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// the total pixel number in the epitome
	int eNum = eHeight*eWidth;
	int pNum = patchSideLengh * patchSideLengh;
	// the total number of the images in the current image category
	int ImageNum = v_imageFileNameList.size();
	// for imread
	int imgReadFlag = 0; // for gray images
	std::string s_numIteration = std::to_string(static_cast<long long>(numIteration));
	std::string fileStorageName = EpitomeResultDir + "/" + whatkindImgs + "_RowCol_Idx.yml";
	// the parameter is convenient to debug 
	unsigned int fftw_flag = FFTW_MEASURE;
	//	unsigned int fftw_flag = FFTW_ESTIMATE;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// some dimension parameters /////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





#define My_Original_FFTW_Convolution
#ifdef  My_Original_FFTW_Convolution

	// ******************************************************
	// ******** My Original FFTW Convolution ****************
	// ******************************************************
	// get the FFTSize, which is equivalent to the size of the convolution of patches (like, 8 by 8) and the Epitome size (like, 50 by 50).
	int FFTSizeWidth = eWidth + patchSideLengh - 1, FFTSizeHeight = eHeight + patchSideLengh - 1;
	// the 2-D input signal length for FFW operation.
	int double_FFTSize_num = FFTSizeWidth * FFTSizeHeight;
	// the output signal length after real-to-complex FFTW.
	int complex_output_num = FFTSizeHeight *(FFTSizeWidth / 2 + 1);


	// before the EM iteration, make a directory for the reconstructed images
	// to define most of the important vector data and fftw_plan before the for-loop
	// often be used vectors by many fftw_plan objects


	typedef vector<double>::iterator v_Double_Iterator;
	typedef vector<int>::iterator v_Int_Iterator;


	std::vector<double>
		v_patchSize_Ones(patchSideLengh*patchSideLengh, 1), // ones
		v_enlarged_Ones(double_FFTSize_num, 0); // because Zero-padding is needed when  v_patchSizeOnes is enlarged to get v_enlarged_Ones, first we set all the values be zeros.
	std::vector<std::complex<double>>
		v_enlarged_OnesFFT(complex_output_num);
	// ones and onesFFT are constant, so they should be calculated before the EM loop
	// to get the v_enlarged_Ones from v_patchSize_Ones, for the following convolution via FFT.

	fftw_plan plan_enlarged_OnesFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_Ones[0], reinterpret_cast<fftw_complex*>(&v_enlarged_OnesFFT[0]), fftw_flag);
	for (int r = 0; r < patchSideLengh; r++){
		for (int c = 0; c < patchSideLengh; c++){
			int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.
			int temp_counter2 = patchSideLengh*r + c;
			// using [] operator
			v_enlarged_Ones[temp_counter] = v_patchSize_Ones[temp_counter2];
		}
	}

	// ones FFT
	fftw_execute(plan_enlarged_OnesFFT);
	// destroy ones FFT plan, which will not delete the onesFFT results in v_enlarged_OnesFFT
	fftw_destroy_plan(plan_enlarged_OnesFFT);
	//delete fft input data
	std::vector<double>().swap(v_patchSize_Ones);
	std::vector<double>().swap(v_enlarged_Ones);
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout << "We are going to learn " << whatkindImgs << " category's Epitome."
		<< "\nThis category has " << ImageNum << " images in total.\n";

	// before the EM algorithm, keep the time
	t = clock() - t;
	printf("It took me %f seconds for initialization before EM.\n", ((float)t) / CLOCKS_PER_SEC);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// EM  Algorithm /////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	long TotalNum = 0;
	for (int em_iterator = 0; em_iterator < numIteration; em_iterator++){ // Beginning of EM Algorithm Iteration

		if (verbose){
			std::cout << "For category: " << whatkindImgs << " ," << em_iterator + 1 << "-th EM Iteration starts now!\n";
		}

		// here we define most of the vector type data and FFTW plans before the EM loop, and let along the xPatch loop.
		std::vector<double>
			v_enlarged_eMeanOverVar(double_FFTSize_num), // will be used during the inner loop, i.e., for each training patch
			v_enlarged_InvVar(double_FFTSize_num), // will be used during the inner loop, i.e., for each training patch
			v_enlarged_LogVar(double_FFTSize_num),
			v_enlarged_eMean2OverVar(double_FFTSize_num);

		std::vector<double> v_double_eLogVarSum(eNum);
		std::vector<double> v_double_eeSum(eNum);

		std::vector<std::complex<double>>
			v_enlarged_eMeanOverVarFFT(complex_output_num), // will be used during the inner loop, i.e., for each training patch
			v_enlarged_LogVarFFT(complex_output_num),
			v_enlarged_TempFFT(complex_output_num),
			v_enlarged_eMean2OverVarFFT(complex_output_num),
			v_enlarged_InvVarFFT(complex_output_num); // will be used during the inner loop, i.e., for each training patch


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// The definition of all the vector data and FFTW plans /////////////
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		// FFT via FFTW
		// since I have put all the elements of the 2-D array in the vector<...>, which makes 1-D FFT available now.
		// If your program performs many transforms of the same size and initialization time is not important, use FFTW_MEASURE; otherwise use FFTW_ESTIMATE.

		fftw_plan
			// for enlarged_LogVarFFT, fft + ifft
			plan_enlarged_LogVarFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_LogVar[0], reinterpret_cast<fftw_complex*>(&v_enlarged_LogVarFFT[0]), fftw_flag),
			i_plan_enlarged_LogVarFFT = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_TempFFT[0]), &v_enlarged_LogVar[0], fftw_flag),

			// for enlarged_eMeanOverVarFFT, only fft
			plan_enlarged_eMeanOverVarFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_eMeanOverVar[0], reinterpret_cast<fftw_complex*>(&v_enlarged_eMeanOverVarFFT[0]), fftw_flag),

			// for enlarged_eMean2OverVarFFT, fft + ifft
			plan_enlarged_eMean2OverVarFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_eMean2OverVar[0], reinterpret_cast<fftw_complex*>(&v_enlarged_eMean2OverVarFFT[0]), fftw_flag),
			i_plan_enlarged_eMean2OverVarFFT = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_TempFFT[0]), &v_enlarged_eMean2OverVar[0], fftw_flag),

			// fft
			plan_enlarged_InvVarFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_InvVar[0], reinterpret_cast<fftw_complex*>(&v_enlarged_InvVarFFT[0]), fftw_flag);


		// initializing the enlarged eMean and eVar, etc.
		// considering the vector data has continuous memory address, thus I will finish their initialization respectively, 
		// that is, not do the initialization just in one for loop.
		// for v_enlarged_LogVar
		for (int r = 0; r < FFTSizeHeight; r++){
			for (int c = 0; c < FFTSizeWidth; c++){
				int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.				
				int tempIdx = (r % eHeight)*eWidth + (c % eWidth);
				double double_enlarged_eMean = v_eMean_Red[tempIdx];
				double double_enlarged_eVar = v_eVar_Red[tempIdx];
				v_enlarged_LogVar[temp_counter] = log(double_enlarged_eVar);
			}
		}

		// for v_enlarged_InvVar
		for (int r = 0; r < FFTSizeHeight; r++){
			for (int c = 0; c < FFTSizeWidth; c++){
				int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.				
				int tempIdx = (r % eHeight)*eWidth + (c % eWidth);
				double double_enlarged_eMean = v_eMean_Red[tempIdx];
				double double_enlarged_eVar = v_eVar_Red[tempIdx];
				v_enlarged_InvVar[temp_counter] = 1 / double_enlarged_eVar;
			}
		}

		// for v_enlarged_eMeanOverVar
		for (int r = 0; r < FFTSizeHeight; r++){
			for (int c = 0; c < FFTSizeWidth; c++){
				int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.				
				int tempIdx = (r % eHeight)*eWidth + (c % eWidth);
				double double_enlarged_eMean = v_eMean_Red[tempIdx];
				double double_enlarged_eVar = v_eVar_Red[tempIdx];
				v_enlarged_eMeanOverVar[temp_counter] = double_enlarged_eMean / double_enlarged_eVar;
			}
		}

		// for v_enlarged_eMean2OverVar
		for (int r = 0; r < FFTSizeHeight; r++){
			for (int c = 0; c < FFTSizeWidth; c++){
				int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.				
				int tempIdx = (r % eHeight)*eWidth + (c % eWidth);
				double double_enlarged_eMean = v_eMean_Red[tempIdx];
				double double_enlarged_eVar = v_eVar_Red[tempIdx];
				v_enlarged_eMean2OverVar[temp_counter] = double_enlarged_eMean * double_enlarged_eMean / double_enlarged_eVar;
			}
		}

		///////////////////////////////////////////////////
		// to calculate v_eLogVarSum /////////////
		///////////////////////////////////////////////////
		fftw_execute(plan_enlarged_LogVarFFT);

		// inverse FFTW plan
		//		plan_enlarged_LogVarFFT = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_TempFFT[0]), &v_enlarged_LogVar[0], fftw_flag);

		// do the multiplication of 2 FFTs
		std::cout << endl;
		for (int temp_counter = 0; temp_counter < complex_output_num; temp_counter++){

			v_enlarged_TempFFT[temp_counter] = v_enlarged_LogVarFFT[temp_counter] * v_enlarged_OnesFFT[temp_counter];
		}


		// inverse FFTW computes an unnormalized DFT. Thus, computing a forward followed by a backward transform(or vice versa) results in the original array scaled by N.
		// if you want to get the original data, keep in mind that the inverse result should be divided by the number N.
		fftw_execute(i_plan_enlarged_LogVarFFT);


		// do the subscripts control of the convolution results
		//		std::vector<double> v_double_eLogVarSum(eNum);
		int temp_counter2 = 0;
		for (int r = patchSideLengh - 1; r < FFTSizeHeight; r++){
			for (int c = patchSideLengh - 1; c < FFTSizeWidth; c++){
				int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.				
				// divided by the length of the input data after doing the inverse FFTW to get the normalized original input signal.
				v_double_eLogVarSum[temp_counter2] = v_enlarged_LogVar[temp_counter] / double_FFTSize_num;
				temp_counter2++;

			}
		}

		// destroy plans which has been executed 
		fftw_destroy_plan(plan_enlarged_LogVarFFT);
		fftw_destroy_plan(i_plan_enlarged_LogVarFFT);
		std::vector<double>().swap(v_enlarged_LogVar);
		std::vector<std::complex<double>>().swap(v_enlarged_LogVarFFT);

		///////////////////////////////////////////////////
		// to calculate v_eeSum /////////////
		///////////////////////////////////////////////////


		fftw_execute(plan_enlarged_eMean2OverVarFFT);

		// inverse FFTW plan
		//		i_plan_enlarged_eMean2OverVarFFT = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_TempFFT[0]), &v_enlarged_eMean2OverVar[0], fftw_flag);

		// do the multiplication of eLogVarFFT and onesFFT
		for (int temp_counter = 0; temp_counter < complex_output_num; temp_counter++){
			v_enlarged_TempFFT[temp_counter] = v_enlarged_eMean2OverVarFFT[temp_counter] * (v_enlarged_OnesFFT[temp_counter]);
		}

		fftw_execute(i_plan_enlarged_eMean2OverVarFFT);

		// do the subscripts control of the convolution results
		//		std::vector<double> v_double_eeSum(eNum);
		temp_counter2 = 0;
		for (int r = patchSideLengh - 1; r < FFTSizeHeight; r++){
			for (int c = patchSideLengh - 1; c < FFTSizeWidth; c++){
				int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.
				// divided by the length of the input data after doing the inverse FFTW to get the normalized original input signal.
				v_double_eeSum[temp_counter2] = v_enlarged_eMean2OverVar[temp_counter] / double_FFTSize_num;
				temp_counter2++;
			}
		}

		// destroy plans which has been executed 
		fftw_destroy_plan(i_plan_enlarged_eMean2OverVarFFT);
		fftw_destroy_plan(plan_enlarged_eMean2OverVarFFT);
		std::vector<double>().swap(v_enlarged_eMean2OverVar);
		std::vector<std::complex<double>>().swap(v_enlarged_eMean2OverVarFFT);
		std::vector<std::complex<double>>().swap(v_enlarged_TempFFT);

		// only fft, but the output should be passed into the Patch loop.
		fftw_execute(plan_enlarged_eMeanOverVarFFT);  // first execute the plan of plan_enlarged_eMeanOverVarFFT to get the v_enlarged_eMeanOverVarFFT
		fftw_execute(plan_enlarged_InvVarFFT); // // first execute the plan of plan_plan_enlarged_InvVarFFT to get the v_enlarged_InvVarFFT
		// destroy plans which has been executed 
		fftw_destroy_plan(plan_enlarged_InvVarFFT);
		fftw_destroy_plan(plan_enlarged_eMeanOverVarFFT);
		std::vector<double>().swap(v_enlarged_eMeanOverVar);
		std::vector<double>().swap(v_enlarged_InvVar);

		// important parameters for calculating the Epitome model
		// before the starting of extracted training patches loop, please make sure zero-initialized doubles to v_double_sumQ(eNum), v_double_sumQX(eNum), and v_double_sumQXX(eNum).
		// i.e., clear the matrices used for collecting sufficient statistics.
		std::vector<double> v_double_sumQ(eNum, 0), v_double_sumQX(eNum, 0), v_double_sumQXX(eNum, 0); // zero-initialized doubles


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// for each image of the current image genre, to extract all the possible training patches for learning the Epitome model
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for (int current_image_index = 0; current_image_index != ImageNum; ++current_image_index){
			Mat mat_currentInputImage = imread(DatabaseDir + "/" + v_imageFileNameList[current_image_index], 0); // read a gray level image
			int channels = mat_currentInputImage.channels();
			// size of input image

			int image_w = mat_currentInputImage.cols;
			int image_h = mat_currentInputImage.rows;

			//vector<int> v_size(2);
			//v_size[0] = image_h;
			//v_size[1] = image_w;
			// show error if no image is read;
			if (mat_currentInputImage.empty())
				cerr << " > No input image is being read!" << endl;


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// generating possible patches extracted out of the current input image
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
			int PatchNum = r_size* c_size;
			int Idx_Num = PatchNum + 2; // the first two elements are image_h and image_w;



			/////////////////////////////////////////////////////////////////////////////////////
			// for each training patch to generate several parameters        //////////////////
			////////////////////////////////////////////////////////////////////////////////////


			// generate the row index and column index for extracting the patch from the input image.
			// vector::end(), Returns an iterator referring to the past-the-end element in the vector container.
			// The past-the-end element is the theoretical element that would follow the last element in the vector. 
			// It does not point to any element, and thus shall be dereferenced.


			////////////////////////////////////////////////////////////////////////////////////
			//************************************************************************************
			//		bool RandomShifting = true; // a parameter just used when debugging process
			//************************************************************************************
			////////////////////////////////////////////////////////////////////////////////////

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


			int trainingCounter = 0; // counting the number of training patches.
			// Here we would calculate vector::end() just once. 
			for (int r_it = 0 ; r_it < r_size; ++ r_it) {
				for (int c_it = 0; c_it < c_size; ++ c_it){
					//				* void srand (unsigned int seed);
					//				* Initialize random number generator
					//				* The pseudo-random number generator is initialized using the argument passed as seed.
					//				* For every different seed value used in a call to srand, the pseudo-random number generator
					//				* can be expected to generate a different succession of results in the subsequent calls to rand.

					// for debugging the code
					//	bool debug_r_it = (r_it == r_end - 1);
					int	temp_r = v_rIdx[r_it];
					int	temp_c = v_cIdx[c_it];

					////////////////////////////////////////////////////////////////////////////////////
					//************************************************************************************
					//	display and estimate the remaining time to finish the process of learning an Epitome model
					//************************************************************************************
					////////////////////////////////////////////////////////////////////////////////////
					// 
					if (verbose & (trainingCounter % (PatchNum / 10) == 0)){
						std::cout << "\nThe position of the current xPatch is: " << "temp_row = " << temp_r << "  and temp_column = " << temp_c << endl;
						std::cout << "Iteration " << em_iterator + 1 << ", image " << current_image_index + 1 << "/" << ImageNum << "  " << (double(100 * trainingCounter)) / PatchNum << "% Complete\n";
						//cout << "Time Remaining:  " << << "seconds\n";
					}

					// do FFT via FFTW
					// do the preparation of implementing FFTW
					// Although I have testing the validity of auto zero-padding by FFTW to input signal when doing FFT and IFFT, I want transform 2-D image signal to 1-D row-major signal,
					// and directly input those 1-D signal to the FFT operator, so I have to do the zero-padding during the 1D-2D transformation, to guarantee the logical correctness.			
					// Hence here we do not use fftw_malloc() to allocate the transformed 1-D input array.


					std::vector<double> //Input Signals
						v_enlarged_xPatch(double_FFTSize_num, .0), // zeros initialization
						v_enlarged_xxPatch(double_FFTSize_num, .0), // zeros initialization
						v_xPatch(patchSideLengh*patchSideLengh, .0),
						v_flipped_xPatch(patchSideLengh*patchSideLengh, .0);
					std::vector<double> v_double_unnormalized_post(eNum);// to calculate the unnormalized posterior		
					std::vector<double> v_double_xxSum(eNum), v_double_xeSum(eNum);

					std::vector<double>
						v_normalized_post(eNum), v_normalized_flipped_post(eNum),
						v_enlarged_flipped_post(double_FFTSize_num, .0), // zero initialization due to the subsequent zero-padding of FFT.
						v_enlarged_Temp(double_FFTSize_num), // to circularly wrap the posterior 
						v_enlarged_TempFlip(double_FFTSize_num); // to hold the values of the flipped version of  v_enlarged_Temp.

					std::vector<std::complex<double>> //FFTs
						v_enlarged_Temp2FFT(complex_output_num),
						v_enlarged_xPatchFFT(complex_output_num),
						v_enlarged_xxPatchFFT(complex_output_num);
					// FFT of Flipped Posterior of mapping m in the Epitome Model.
					std::vector<complex<double>> v_Flipped_PostFFT(complex_output_num);

					// fft + ifft
					fftw_plan
						plan_enlarged_xPatchFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_xPatch[0], reinterpret_cast<fftw_complex*>(&v_enlarged_xPatchFFT[0]), fftw_flag),
						plan_enlarged_xxPatchFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_xxPatch[0], reinterpret_cast<fftw_complex*>(&v_enlarged_xxPatchFFT[0]), fftw_flag),
						i_plan_enlarged_xPatchFFT = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_Temp2FFT[0]), &v_enlarged_xPatch[0], fftw_flag),
						i_plan_enlarged_xxPatchFFT = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_Temp2FFT[0]), &v_enlarged_xxPatch[0], fftw_flag),

						// fft + ifft
						plan_sumQ = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_flipped_post[0], reinterpret_cast<fftw_complex*>(&v_Flipped_PostFFT[0]), fftw_flag),
						i_plan_sumQ = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_Temp2FFT[0]), &v_enlarged_Temp[0], fftw_flag),

						// only ifft
						i_plan_sumQX = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_Temp2FFT[0]), &v_enlarged_Temp[0], fftw_flag),
						i_plan_sumQXX = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_Temp2FFT[0]), &v_enlarged_Temp[0], fftw_flag);

					///////////////////////////////////////////////////////////////////////////////
					/////////// for each possible xPatchFFT and xxPatchFFT, doing E-Step ////////////////////////
					///////////////////////////////////////////////////////////////////////////////


					// casting STL complex<double> to fftw_complex for FFT and IFFT via the FFTW library.
					// Syntax: std::vector<std::complex<double> > a1(N), a2(N); // N complex values in the input and output data
					// fftw_plan_dft(N, reinterpret_cast<fftw_complex*>(&a1[0]),reinterpret_cast<fftw_complex*>(&a2[0]), FFTW_FORWARD, FFTW_ESTIMATE);	




					// to get the extracted Patch for the following convolution via FFT.
					temp_counter2 = 0;
					for (int i = temp_r; i < temp_r + patchSideLengh; i++) {
						for (int j = temp_c; j < temp_c + patchSideLengh; j++){
							// initializing the element of the 2-D array.
							if (channels == 1){// gray-level image
								if ((i > image_h) | (j > image_w)){
									cerr << "Indices Exceeding is happening!" << endl;
								}
								v_xPatch[temp_counter2] = (double)(mat_currentInputImage.at<uchar>(i, j)) / 255.0f;
								temp_counter2++;
							}
						}
					}


					// to get the flipped xPatch for the following convolution via FFT.

					temp_counter2 = 0;
					for (int r = patchSideLengh - 1; r >= 0; r--){
						for (int c = patchSideLengh - 1; c >= 0; c--){
							int temp_counter = patchSideLengh*r + c; //2-D array indexing for the element in r-th row and c-th column.
							v_flipped_xPatch[temp_counter2] = v_xPatch[temp_counter];
							temp_counter2++;
						}
					}



					// r2c FFTW plan
					//					plan_enlarged_xPatchFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_xPatch[0], reinterpret_cast<fftw_complex*>(&v_enlarged_xPatchFFT[0]), fftw_flag);
					//					plan_enlarged_xxPatchFFT = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_xxPatch[0], reinterpret_cast<fftw_complex*>(&v_enlarged_xxPatchFFT[0]), fftw_flag);

					// to get the enlarged_xPatch from the flipped xPatch, for the following convolution via FFT.
					for (int r = 0; r < patchSideLengh; r++){
						for (int c = 0; c < patchSideLengh; c++){
							int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.
							int temp_counter2 = patchSideLengh*r + c;

							v_enlarged_xPatch[temp_counter] = v_flipped_xPatch[temp_counter2];

						}
					}

					// to get the enlarged_xxPatch from the flipped xPatch, for the following convolution via FFT.
					for (int r = 0; r < patchSideLengh; r++){
						for (int c = 0; c < patchSideLengh; c++){
							int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.
							int temp_counter2 = patchSideLengh*r + c;

							v_enlarged_xxPatch[temp_counter] = v_flipped_xPatch[temp_counter2] * v_flipped_xPatch[temp_counter2];

						}
					}
					// doing forward xPatch and xxPatch FFT
					fftw_execute(plan_enlarged_xPatchFFT);
					fftw_execute(plan_enlarged_xxPatchFFT);


					// build inverse FFTW plans before doing multiplication of two FFTs
					// You must create the plan before initializing the input, because the FFTW_MEASURE overwrites the in/out arrays.
					// Initializing the FFTW input, here the number of the element which I will initialize is patchSideLengh*patchSideLengh, and the resting elements will be initialized as zeros, 
					// due to the auto zero-padding operation of FFTW.

					// IFFTW
					//					i_plan_enlarged_xPatchFFT = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_Temp2FFT[0]), &v_enlarged_xPatch[0], fftw_flag);
					//					i_plan_enlarged_xxPatchFFT = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_Temp2FFT[0]), &v_enlarged_xxPatch[0], fftw_flag);


					for (int temp_counter = 0; temp_counter < complex_output_num; temp_counter++){

						v_enlarged_Temp2FFT[temp_counter] = v_enlarged_eMeanOverVarFFT[temp_counter] * v_enlarged_xPatchFFT[temp_counter];
					}

					// doing inverse FFT via FFTW
					fftw_execute(i_plan_enlarged_xPatchFFT);


					for (int temp_counter = 0; temp_counter < complex_output_num; temp_counter++){

						v_enlarged_Temp2FFT[temp_counter] = v_enlarged_InvVarFFT[temp_counter] * v_enlarged_xxPatchFFT[temp_counter];

					}

					// doing inverse FFT via FFTW
					fftw_execute(i_plan_enlarged_xxPatchFFT);


					///////////////////////////////////////////////////////////////////////////////
					/////////// for xxSum  //////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////

					// do the subscripts control of the convolution results
					//					std::vector<double> v_double_xxSum(eNum), v_double_xeSum(eNum);
					temp_counter2 = 0;
					for (int r = patchSideLengh - 1; r < FFTSizeHeight; r++){
						for (int c = patchSideLengh - 1; c < FFTSizeWidth; c++){
							int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.		
							// doing normalization after inverse FFT via FFTW				
							v_double_xxSum[temp_counter2] = v_enlarged_xxPatch[temp_counter] / double_FFTSize_num;
							temp_counter2++;
						}
					}

					///////////////////////////////////////////////////////////////////////////////
					/////////// for xeSum  //////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////

					// do the subscripts control of the convolution results
					//					std::vector<double> v_double_xxSum(eNum), v_double_xeSum(eNum);
					temp_counter2 = 0;
					for (int r = patchSideLengh - 1; r < FFTSizeHeight; r++){
						for (int c = patchSideLengh - 1; c < FFTSizeWidth; c++){
							// static int temp_counter2 = 0;
							int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.

							// doing normalization after inverse FFT via FFTW
							v_double_xeSum[temp_counter2] = v_enlarged_xPatch[temp_counter] / double_FFTSize_num;

							temp_counter2++;
						}
					}

					// destroy fft plans and their unused corresponding vector data
					fftw_destroy_plan(i_plan_enlarged_xxPatchFFT);
					fftw_destroy_plan(i_plan_enlarged_xPatchFFT);
					fftw_destroy_plan(plan_enlarged_xxPatchFFT);
					fftw_destroy_plan(plan_enlarged_xPatchFFT);
					std::vector<double>().swap(v_xPatch);
					std::vector<double>().swap(v_flipped_xPatch);


					// to calculate the unnormalized posterior or likelihood matrix,
					// and find the maximum element
					double max_unnormalized_post = DOUBLE_MIN; // to get a minimum value in double
					for (int r = 0; r < eHeight; r++){
						for (int c = 0; c < eWidth; c++){
							int temp_counter = eWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.
							// using operator [] to access element of vector
							v_double_unnormalized_post[temp_counter] = (-0.5)* (v_double_eeSum[temp_counter] - 2 * v_double_xeSum[temp_counter] + v_double_xxSum[temp_counter] + v_double_eLogVarSum[temp_counter]);

							if (max_unnormalized_post < v_double_unnormalized_post[temp_counter]){
								max_unnormalized_post = v_double_unnormalized_post[temp_counter];
							}
						}
					}



					///////////////////////////////////////////////////////////////////////////////
					/////////// for each possible xPatchFFT and xxPatchFFT, doing M-Step ////////////////////////
					///////////////////////////////////////////////////////////////////////////////

					// to normalized the above posterior
					// to introduce a constant (named alpha) to avoid the overflow error during do exp(...)
					double alpha = max_unnormalized_post - 0.5*log(DOUBLE_MAX) + 2 * log(double (patchSideLengh ^ 2));
					double sum_exp_post = 0;
					for (int temp_counter = 0; temp_counter < eNum; temp_counter++){
						// using operator [] to access element of vector
						sum_exp_post += exp(v_double_unnormalized_post[temp_counter] - alpha);
					}


					double normalizingConst = log(sum_exp_post) + alpha;
					/*
					std::vector<double>
					v_normalized_post(eNum), v_normalized_flipped_post(eNum),
					v_enlarged_flipped_post(double_FFTSize_num, 0), // zero initialization due to the subsequent zero-padding of FFT.
					v_enlarged_Temp(double_FFTSize_num), // to circularly wrap the posterior
					v_enlarged_TempFlip(double_FFTSize_num); // to hold the values of the flipped version of  v_enlarged_Temp.
					*/
					// normalized posterior
					for (int temp_counter = 0; temp_counter < eNum; temp_counter++){
						// using operator [] to access element of vector
						v_normalized_post[temp_counter] = exp(v_double_unnormalized_post[temp_counter] - normalizingConst);
					}


					// to get the flipped normalized posterior, since xPatch and xxPatch have already been flipped before.
					temp_counter2 = 0;
					for (int r = eHeight - 1; r >= 0; r--){
						for (int c = eWidth - 1; c >= 0; c--){
							int temp_counter = eWidth *r + c; //2-D array indexing for the element in r-th row and c-th column.
							// using operator [] to access element of vector
							v_normalized_flipped_post[temp_counter2] = v_normalized_post[temp_counter];
							temp_counter2++;
						}
					}


					// FFT of Flipped Posterior of mapping m in the Epitome Model.
					//					std::vector<complex<double>> v_Flipped_PostFFT(complex_output_num);

					// r2c FFTW plan
					//					plan_sumQ = fftw_plan_dft_r2c_1d(double_FFTSize_num, &v_enlarged_flipped_post[0], reinterpret_cast<fftw_complex*>(&v_Flipped_PostFFT[0]), fftw_flag);

					// initialize the enlarged flipped posteriors
					for (int r = 0; r < eHeight; r++){
						for (int c = 0; c < eWidth; c++){
							int temp_counter = FFTSizeWidth*r + c; //2-D array indexing for the element in r-th row and c-th column.
							int temp_counter2 = eWidth*r + c;
							// using operator [] to access element of vector
							v_enlarged_flipped_post[temp_counter] = v_normalized_flipped_post[temp_counter2];
						}
					}



					///////////////////////////////////////////////////////////////////////////////
					/////////// to calculate v_double_sumQ, v_double_sumQX,v_double_sumQXX ////////////////////////
					///////////////////////////////////////////////////////////////////////////////

					/////////// to calculate v_double_sumQ ////////////

					// Do not flip the posterior when computing the FFT, since this time we need convolutions instead of correlations.
					// But pay atention to use the non-flipped xPatch and xxpatch;
					// or if you want to use the flipped xPatch and xxPatch, please flip or even re-flip the posterior appropriately.

					// do FFT for calculating sumQ
					fftw_execute(plan_sumQ);
					// IFFT plan
					//					i_plan_sumQ = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_TempFFT[0]), &v_enlarged_Temp[0], fftw_flag);

					// do multiplication of 2 FFTs
					for (int temp = 0; temp < complex_output_num; temp++){
						v_enlarged_Temp2FFT[temp] = v_Flipped_PostFFT[temp] * v_enlarged_OnesFFT[temp];
					}
					// do IFFT
					fftw_execute(i_plan_sumQ);

					// to flip v_enlarged_Temp
					temp_counter2 = 0;
					for (int r = FFTSizeHeight - 1; r >= 0; r--){
						for (int c = FFTSizeWidth - 1; c >= 0; c--){
							int temp_counter = FFTSizeWidth *r + c; //2-D array indexing for the element in r-th row and c-th column.
							v_enlarged_TempFlip[temp_counter2] = v_enlarged_Temp[temp_counter];
							temp_counter2++;
						}
					}

					// circularly wrap, for two column blocks 
					// e.g., patchSideLength = 8, we select the the first 7-column block and the last 7-column block, and send the sum of the two blocks to the first 7-column block.
					for (int r = 0; r < FFTSizeHeight; r++){
						for (int c = 0; c < patchSideLengh - 1; c++){
							int temp_counter = FFTSizeWidth*r + c; //first 7-column block.
							int temp_counter2 = FFTSizeWidth*r + (eWidth + c); //the last 7-column block
							// using operator [] to access element of vector
							v_enlarged_TempFlip[temp_counter] += v_enlarged_TempFlip[temp_counter2];
						}
					}

					// circularly wrap, for two row blocks 
					// E.G. patchSideLength = 8, we select the the first 7-row block and the last 7-row block, and send the sum of the two row blocks to the first 7-row block.
					for (int r = 0; r < patchSideLengh - 1; r++){
						for (int c = 0; c < FFTSizeWidth; c++){
							int temp_counter = FFTSizeWidth*r + c; //first 7-row block.
							int temp_counter2 = FFTSizeWidth*(r + eWidth) + c; //the last 7-row block
							// using operator [] to access element of vector
							v_enlarged_TempFlip[temp_counter] += v_enlarged_TempFlip[temp_counter2];
						}
					}

					// do the subscripts control of the convolution results
					for (int r = 0; r < eHeight; r++){
						for (int c = 0; c < eWidth; c++){
							int temp = r* FFTSizeWidth + c;
							int temp2 = r* eWidth + c;
							// using operator [] to access element of vector
							// "+=" is needed, instead of "+", because "+=" can guarantee to get the sum of parameter calculation of each training patch.
							// normalization after the inverse FFT via the FFTW
							v_double_sumQ[temp2] += v_enlarged_TempFlip[temp] / double_FFTSize_num;
						}
					}


					/////////// to calculate v_double_sumQX ////////////

					// IFFT, sumQX FFTW plan
					//					i_plan_sumQX = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_TempFFT[0]), &v_enlarged_Temp[0], fftw_flag);


					// do multiplication of 2 FFTs
					for (int temp = 0; temp < complex_output_num; temp++){
						// using operator [] to access element of vector
						v_enlarged_Temp2FFT[temp] = v_Flipped_PostFFT[temp] * v_enlarged_xPatchFFT[temp];
					}

					// do IFFT
					fftw_execute(i_plan_sumQX);

					// flip the convolution result, i.e., to flip v_enlarged_Temp
					temp_counter2 = 0;
					for (int r = FFTSizeHeight - 1; r >= 0; r--){
						for (int c = FFTSizeWidth - 1; c >= 0; c--){
							int temp_counter = FFTSizeWidth *r + c; //2-D array indexing for the element in r-th row and c-th column.
							// using operator [] to access element of vector
							v_enlarged_TempFlip[temp_counter2] = v_enlarged_Temp[temp_counter];
							temp_counter2++;
						}
					}

					// circularly wrap, for two column blocks 
					// e.g., patchSideLength = 8, we select the the first 7-column block and the last 7-column block, and send the sum of the two blocks to the first 7-column block.
					for (int r = 0; r < FFTSizeHeight; r++){
						for (int c = 0; c < patchSideLengh - 1; c++){
							int temp_counter = FFTSizeWidth*r + c; //first 7-column block.
							int temp_counter2 = FFTSizeWidth*r + (eWidth + c); //the last 7-column block
							// using operator [] to access element of vector
							v_enlarged_TempFlip[temp_counter] += v_enlarged_TempFlip[temp_counter2];
						}
					}

					// circularly wrap, for two row blocks 
					// e.g., patchSideLength = 8, we select the the first 7-row block and the last 7-row block, and send the sum of the two row blocks to the first 7-row block.
					for (int r = 0; r < patchSideLengh - 1; r++){
						for (int c = 0; c < FFTSizeWidth; c++){
							int temp_counter = FFTSizeWidth*r + c; //first 7-row block.
							int temp_counter2 = FFTSizeWidth*(r + eWidth) + c; //the last 7-row block
							// using operator [] to access element of vector
							v_enlarged_TempFlip[temp_counter] += v_enlarged_TempFlip[temp_counter2];
						}
					}

					// do the subscripts control of the convolution results
					for (int r = 0; r < eHeight; r++){
						for (int c = 0; c < eWidth; c++){
							int temp = r* FFTSizeWidth + c;
							int temp2 = r* eWidth + c;
							// using operator [] to access element of vector
							// "+=" is needed, instead of "+", because "+=" can guarantee to get the sum of parameter calculation of each training patch.
							// normalization after the inverse FFT via the FFTW
							v_double_sumQX[temp2] += v_enlarged_TempFlip[temp] / double_FFTSize_num;
						}
					}


					/////////// to calculate v_double_sumQXX ////////////
					// IFFT, sumQXX FFTW plan
					//					i_plan_sumQXX = fftw_plan_dft_c2r_1d(double_FFTSize_num, reinterpret_cast<fftw_complex*>(&v_enlarged_TempFFT[0]), &v_enlarged_Temp[0], fftw_flag);

					// do multiplication of 2 FFTs
					for (int temp = 0; temp < complex_output_num; temp++){
						// using operator [] to access element of vector
						v_enlarged_Temp2FFT[temp] = v_Flipped_PostFFT[temp] * v_enlarged_xxPatchFFT[temp];
					}

					// do IFFT
					fftw_execute(i_plan_sumQXX);

					// flip the convolution result
					// to flip v_enlarged_Temp
					temp_counter2 = 0;
					for (int r = FFTSizeHeight - 1; r >= 0; r--){
						for (int c = FFTSizeWidth - 1; c >= 0; c--){
							int temp_counter = FFTSizeWidth *r + c; //2-D array indexing for the element in r-th row and c-th column.
							// using operator [] to access element of vector
							v_enlarged_TempFlip[temp_counter2] = v_enlarged_Temp[temp_counter];
							temp_counter2++;
						}
					}

					// circularly wrap, for two column blocks 
					// e.g., patchSideLength = 8, we select the the first 7-column block and the last 7-column block, and send the sum of the two blocks to the first 7-column block.
					for (int r = 0; r < FFTSizeHeight; r++){
						for (int c = 0; c < patchSideLengh - 1; c++){
							int temp_counter = FFTSizeWidth*r + c; //first 7-column block.
							int temp_counter2 = FFTSizeWidth*r + (eWidth + c); //the last 7-column block
							// using operator [] to access element of vector
							v_enlarged_TempFlip[temp_counter] += v_enlarged_TempFlip[temp_counter2];
						}
					}

					// circularly wrap, for two row blocks 
					// e.g., patchSideLength = 8, we select the the first 7-row block and the last 7-row block, and send the sum of the two row blocks to the first 7-row block.
					for (int r = 0; r < patchSideLengh - 1; r++){
						for (int c = 0; c < FFTSizeWidth; c++){
							int temp_counter = FFTSizeWidth*r + c; //first 7-row block.
							int temp_counter2 = FFTSizeWidth*(r + eWidth) + c; //the last 7-row block
							// using operator [] to access element of vector
							v_enlarged_TempFlip[temp_counter] += v_enlarged_TempFlip[temp_counter2];
						}
					}

					// do the subscripts control of the convolution results
					for (int r = 0; r < eHeight; r++){
						for (int c = 0; c < eWidth; c++){
							int temp = r* FFTSizeWidth + c;
							int temp2 = r* eWidth + c;
							// using operator [] to access element of vector
							// "+=" is needed, instead of "+", because "+=" can guarantee to get the sum of parameter calculation of each training patch.
							// normalization after the inverse FFT via the FFTW
							v_double_sumQXX[temp2] += v_enlarged_TempFlip[temp] / double_FFTSize_num;
						}
					}

					trainingCounter += 1; // // counting the number of training patches.
					TotalNum++;

					//  delete FFTW plan objects, and to release the memory
					fftw_destroy_plan(i_plan_sumQXX);
					fftw_destroy_plan(i_plan_sumQX);
					fftw_destroy_plan(i_plan_sumQ);
					fftw_destroy_plan(plan_sumQ);
					// delete vector data, free memory
					std::vector<double>().swap(v_enlarged_xPatch);
					std::vector<double>().swap(v_enlarged_xxPatch);
					std::vector<double>().swap(v_double_unnormalized_post);
					std::vector<double>().swap(v_double_xxSum);
					std::vector<double>().swap(v_double_xeSum);
					std::vector<double>().swap(v_normalized_post);
					std::vector<double>().swap(v_normalized_flipped_post);
					std::vector<double>().swap(v_enlarged_flipped_post);
					std::vector<double>().swap(v_enlarged_Temp);
					std::vector<double>().swap(v_enlarged_TempFlip);
					std::vector<std::complex<double>>().swap(v_enlarged_Temp2FFT);
					std::vector<std::complex<double>>().swap(v_enlarged_xPatchFFT);
					std::vector<std::complex<double>>().swap(v_enlarged_xxPatchFFT);
					std::vector<std::complex<double>>().swap(v_Flipped_PostFFT);


				} // End of Patches Extraction and Training (for column indexing)

			}  // End of Patches Extraction and Training (for row indexing)

			if (verbose){
				std::cout << "\n The Current Input Image: " << current_image_index + 1 << "/" << ImageNum << " has been read!\n";
			}


			// delete possible indices of each xPatch in currently read image
			std::vector<int>().swap(v_rIdx);
			std::vector<int>().swap(v_cIdx);

		} // end of each current input image

		// avoid numerical problems
		// to find the nearly -being-zero elements in sumQ
		for (int r = 0; r < eHeight; r++){
			for (int c = 0; c < eWidth; c++){
				int temp = c + r*eWidth;
				// using operator [] to access element of vector
				if (v_double_sumQ[temp] <= tolerance){
					v_double_sumQ[temp] = 1;
					v_double_sumQX[temp] = v_eMean_Red[temp];
					v_double_sumQXX[temp] = v_eVar_Red[temp] + v_eMean_Red[temp] * v_eMean_Red[temp];
				}
			}
		}


		// compute the new epitome
		for (int r = 0; r < eHeight; r++){
			for (int c = 0; c < eWidth; c++){
				int temp = c + r*eWidth;
				// using operator [] to access element of vector
				v_eMean_Red[temp] = v_double_sumQX[temp] / v_double_sumQ[temp];
				v_eVar_Red[temp] = v_double_sumQXX[temp] / v_double_sumQ[temp] - v_eMean_Red[temp] * v_eMean_Red[temp];

				// eMean belongs to [0, 1] interval.
				// make sure that the mean is within the range 0 - 1 (may not because of numerical issues)
				v_eMean_Red[temp] = v_eMean_Red[temp] > 0 ? v_eMean_Red[temp] : 0;
				v_eMean_Red[temp] = v_eMean_Red[temp] < 1 ? v_eMean_Red[temp] : 1;
				//enforce a minimum variance
				v_eVar_Red[temp] = v_eVar_Red[temp] > minVar ? v_eVar_Red[temp] : minVar;

			}
		}

		std::cout << "\nThe " << em_iterator + 1 << "-th EM Iteration is finished now!\n";
		//  delete FFTW plan objects, and to release the memory
		// delete vector data, free memory
		std::vector<std::complex<double>>().swap(v_enlarged_eMeanOverVarFFT);
		std::vector<std::complex<double>>().swap(v_enlarged_InvVarFFT);
		std::vector<double>().swap(v_double_sumQ);
		std::vector<double>().swap(v_double_sumQX);
		std::vector<double>().swap(v_double_sumQXX);
		std::vector<double>().swap(v_double_eLogVarSum);
		std::vector<double>().swap(v_double_eeSum);


	} // End of EM Algorithm


#endif //  My_Original_FFTW_Convolution

	//  delete FFTW plan objects, and to release the memory
	// delete vector data, free memory, which was built before the EM-loop, so it is removed after the finishing of EM loop.
	std::vector<std::complex<double>>().swap(v_enlarged_OnesFFT);

	Mat_<double> mat_eMean(eHeight, eWidth, .0);
	Mat_<double> mat_eVar(eHeight, eWidth, .0);
	// save the vector eMean and eVar into to mat_eMean and mat_eVar for the following txtFile saving
	for (int r = 0; r < eHeight; r++){
		for (int c = 0; c < eWidth; c++){

			// using operator [] to access element of vector
			int temp = r*eWidth + c;
			mat_eMean.at<double>(r, c) = v_eMean_Red[temp];
			mat_eVar.at<double>(r, c) = v_eVar_Red[temp];
		}

	}

	// save mat_eMean into to an image
	cv::imwrite(EpitomeResultDir + "/" + whatkindImgs + "_eMean_EM" + s_numIteration + imgEncodeType, 255 * mat_eMean);

	// save the vector eMean into ".yml" file for image reconstruction
	FileStorage fs_Epitome_in(EpitomeResultDir + "/" + whatkindImgs + "_eMean_EM" + s_numIteration + ".yml", FileStorage::WRITE);
	if (!fs_Epitome_in.isOpened())
	{
		cerr << "failed to open " << EpitomeResultDir + "/" + whatkindImgs + "_eMean_EM" + s_numIteration + ".yml" << endl;
		fileStorageHelp();
	}

	fs_Epitome_in << "e_Mean" << v_eMean_Red;
	// eVaar is not necessary not the following image reconstruction
	//	fs_Epitome_in << "e_Var" << v_eVar;
	fs_Epitome_in.release();

	if (verbose){
		std::cout << " The newly learned Epitome has been saved to mat_eMean and mat_eVar.\n"
			<< " The newly learned mat_eMean has been saved as a jpg image.\n"
			<< " The newly learned mat_eMean and mat_eVar has been saved into " + whatkindImgs + "_eMeanVar_EM" + s_numIteration + ".yml\n"
			<< " EM has finished! TotalNum = " << TotalNum << "\n";
	}

	// to be continued .......
	// to be continued .......
	// to be continued .......

}