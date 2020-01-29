#include "ImageEpitome.hpp"

void ImgEpitome::getMaxPostRowColIdx(
	const vector<string> & v_imageFileNameList, // read input images via their file-names.
	const int & patchSpac){

	// print the consuming time
	clock_t t = clock();
	// for imread
	int imgReadFlag = 0;
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// some dimension parameters /////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	// the total pixel number in the epitome
	int eNum = eHeight*eWidth;
	int pNum = patchSideLengh * patchSideLengh;
	// the total number of the images in the current image category
	int ImageNum = v_imageFileNameList.size();
	
	// to get possible indices of each xPatch in currently read image
	// set new value to the patchSpacing, for example:
	// 1). let patchSpacing = patchSideLengh, i.e., non-overlapping extraction of xPatch;
	// 2). let patchSpacing < patchSideLengh, i.e., overlapping extraction of xPatch;
	this->patchSpacing = patchSpac;
	

	std::string s_numIteration = std::to_string(static_cast<long long>(numIteration));
	std::string fileStorageName = EpitomeResultDir + "/" 
		+ whatkindImgs + "_" + to_string(static_cast<long long>(patchSpacing)) 
		+ "_" + to_string(static_cast<long long>(Read_eMean_via_File_Flag)) + "_RowCol_Idx.yml";

	cout << "We are going to learn " << whatkindImgs << " category's Epitome."
		<< "\nThis category has " << ImageNum << " images in total.\n";

	//*******************************************************
	// just write the beginning of the FileStorage file.
	//*******************************************************
	FileStorage fs_in(fileStorageName, FileStorage::WRITE);
	if (!fs_in.isOpened())
	{
		cerr << "failed to open " << fileStorageName << endl;
		fileStorageHelp();
	}

	fs_in << whatkindImgs << "[:";
	//fs_in.release();

	// before getting the indices, keep the time first
	t = clock() - t;
	printf("It took me %f seconds for initialization before learning the indices.\n", ((float)t) / CLOCKS_PER_SEC);

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

	// Src Signal
	vector<Epitome_Doulbe> v_eMeanOverVar(eNum);
	vector<Epitome_Doulbe> v_InvVar(eNum);
	vector<Epitome_Doulbe> v_LogVar(eNum); // not xPatch
	vector<Epitome_Doulbe> v_eMean2OverVar(eNum); // not xPatch

	// FFT results
	int Num_FFT_Complex_Outputs = ws.h_fftw * (ws.w_fftw / 2 + 1);


	//kernel
	vector<Epitome_Doulbe> v_xPatch(pNum);
	vector<Epitome_Doulbe> v_xxPatch(pNum);
	vector<Epitome_Doulbe> v_patchSize_Ones(pNum, 1); // ones

	// Convolution or Correlation Result
	vector<Epitome_Doulbe> v_eLogVarSum(eNum);
	vector<Epitome_Doulbe> v_eeSum(eNum);
	vector<Epitome_Doulbe> v_double_unnormalized_post(eNum); // to calculate the unnormalized posterior
	vector<Epitome_Doulbe> v_double_xxSum(eNum);
	vector<Epitome_Doulbe> v_double_xeSum(eNum);


	// EM starts here
	int em_iterator = 1;

		// for v_LogVar
		for (int i = 0; i != eNum; ++i)
			v_LogVar[i] = log(v_eVar_Red[i]);

		// for v_InvVar
		for (int i = 0; i != eNum; ++i)
			v_InvVar[i] = 1.0 / v_eVar_Red[i];

		// for v_eMeanOverVar
		for (int i = 0; i != eNum; ++i)
			v_eMeanOverVar[i] = v_eMean_Red[i] / v_eVar_Red[i];

		// for v_eMean2OverVar
		for (int i = 0; i != eNum; ++i)
			v_eMean2OverVar[i] = v_eMean_Red[i] * v_eMean_Red[i] / v_eVar_Red[i];

		///////////////////////////////////////////////////
		// to calculate v_eLogVarSum /////////////
		///////////////////////////////////////////////////
		ws.mode = c_fft::LINEAR_CORRELATION_VALIDSAME;
		update_workspace_few(ws);
		//	c_fft::convolve(ws, v_LogVar, v_patchSize_Ones);
		c_fft::convolve(ws, &v_LogVar[0], &v_patchSize_Ones[0]);
		for (int i = 0; i != eNum; ++i)
			v_eLogVarSum[i] = ws.dst_convolution[i];

		///////////////////////////////////////////////////
		// to calculate v_eeSum /////////////
		///////////////////////////////////////////////////
		//	c_fft::convolve(ws, v_eMean2OverVar, v_patchSize_Ones);
		c_fft::convolve(ws, &v_eMean2OverVar[0], &v_patchSize_Ones[0]);
		for (int i = 0; i != eNum; ++i)
			v_eeSum[i] = ws.dst_convolution[i];


		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// for each image of the current image genre, to extract all the training patches for reconstruction via the already learned Epitome model
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		for (int current_image_index = 0; current_image_index != ImageNum; ++current_image_index){
			Mat mat_currentInputImage = imread(DatabaseDir + "/" + v_imageFileNameList[current_image_index], imgReadFlag); // read a gray level image
			int channels = mat_currentInputImage.channels();
			// size of input image

			int image_w = mat_currentInputImage.cols;
			int image_h = mat_currentInputImage.rows;

			vector<int> v_size(2);
			v_size[0] = image_h;
			v_size[1] = image_w;
			// show error if no image is read;
			if (mat_currentInputImage.empty())
				cerr << " > No input image is being read!" << endl;


			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			// generating possible patches extracted out of the current input image
			// last possible location of patches in the input image x;
			int c_xEndIdx = image_w - patchSideLengh,
				r_xEndIdx = image_h - patchSideLengh;

			
			tuple<vector<int>, vector<int>> t_v_row_col_Idx = getRowColIdxOfOneImg(image_h, image_w);
			vector<int> v_rIdx = get<0>(t_v_row_col_Idx);
			vector<int> v_cIdx = get<1>(t_v_row_col_Idx);

			int r_size = v_rIdx.size(), c_size = v_cIdx.size();
			int PatchNum = r_size* c_size;
			int Idx_Num = PatchNum + 2; // the first two elements are image_h and image_w;

			// the corresponding row_column indices of the maximum posterior element for all the possible xPatch extracted out of the current image
				v_maxPost_row_column_Idx[current_image_index] = vector<int>(Idx_Num, 0);
				v_maxPost_row_column_Idx[current_image_index][0] = image_h;
				v_maxPost_row_column_Idx[current_image_index][1] = image_w;

			/////////////////////////////////////////////////////////////////////////////////////
			// for each training patch to generate several parameters        //////////////////
			////////////////////////////////////////////////////////////////////////////////////
            
			// to reconstruct images, no random shifting and wiggling are needed.
			int trainingCounter = 0; // counting the number of training patches.
			for (int r_it = 0; r_it < r_size; ++r_it) {
				for (int c_it = 0; c_it < c_size; ++c_it){
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
						std::cout <<"image " << current_image_index + 1 << "/" << ImageNum << "  " << (double(100 * trainingCounter)) / PatchNum << "% Complete\n";
					}

					///////////////////////////////////////////////////////////////////////////////
					/////////// for each possible xPatchFFT and xxPatchFFT, doing E-Step ////////////////////////
					///////////////////////////////////////////////////////////////////////////////


					// casting STL complex<double> to fftw_complex for FFT and IFFT via the FFTW library.
					// Syntax: std::vector<std::complex<double> > a1(N), a2(N); // N complex values in the input and output data
					// fftw_plan_dft(N, reinterpret_cast<fftw_complex*>(&a1[0]),reinterpret_cast<fftw_complex*>(&a2[0]), FFTW_FORWARD, FFTW_ESTIMATE);	



					// to get the extracted Patch v_xPatch.
					int temp_counter = 0;
					for (int i = temp_r; i < temp_r + patchSideLengh; ++i)
						for (int j = temp_c; j < temp_c + patchSideLengh; ++j, ++temp_counter)
							v_xPatch[temp_counter] = (double)(mat_currentInputImage.at<uchar>(i, j)) / 255.0f;

					// to get the extracted Patch v_xxPatch.
					temp_counter = 0;
					for (int i = temp_r; i < temp_r + patchSideLengh; i++)
						for (int j = temp_c; j < temp_c + patchSideLengh; j++, ++temp_counter)
							v_xxPatch[temp_counter] = v_xPatch[temp_counter] * v_xPatch[temp_counter];


					///////////////////////////////////////////////////////////////////////////////
					/////////// for xxSum  //////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////
					// correlation
					ws.mode = c_fft::LINEAR_CORRELATION_VALIDSAME;
					update_workspace_few(ws);
					c_fft::convolve(ws, &v_InvVar[0], &v_xxPatch[0]);
					for (int i = 0; i != eNum; ++i)
						v_double_xxSum[i] = ws.dst_convolution[i];

					///////////////////////////////////////////////////////////////////////////////
					/////////// for xeSum  //////////////////////////////////////////////
					///////////////////////////////////////////////////////////////////////////////
					//			c_fft::convolve(ws, v_eMeanOverVar, v_xPatch);
					c_fft::convolve(ws, &v_eMeanOverVar[0], &v_xPatch[0]);
					for (int i = 0; i != eNum; ++i)
						v_double_xeSum[i] = ws.dst_convolution[i];


					// to calculate the unnormalized posterior or likelihood matrix,
					// and find the maximum element wherein and its corresponding index

					double max_unnormalized_post = DOUBLE_MIN; // to get a minimum value in double
					int r_maxPosterior = 0, c_maxPosterior = 0;
					for (int r = 0; r != eHeight; ++r){
						for (int c = 0; c != eWidth; ++c){
							int i = r *eWidth + c;
							v_double_unnormalized_post[i] = (-0.5)* (v_eeSum[i] - 2 * v_double_xeSum[i] + v_double_xxSum[i] + v_eLogVarSum[i]);
							if (max_unnormalized_post < v_double_unnormalized_post[i]){
								max_unnormalized_post = v_double_unnormalized_post[i];
								r_maxPosterior = r;
								c_maxPosterior = c;
							}
						}
					}

					//////////////////////////////////////////////////////////////////////////////////
					// save the corresponding row_column indices of the maximum posterior element 
					// for all the possible xPatch extracted out of the current image
					//////////////////////////////////////////////////////////////////////////////////
					v_maxPost_row_column_Idx[current_image_index][trainingCounter + 2] = r_maxPosterior* eWidth + c_maxPosterior;
					
					trainingCounter += 1; // // counting the number of training patches.

				} // End of Patches Extraction and Training (for column indexing)

			}  // End of Patches Extraction and Training (for row indexing)

			if (verbose){
				std::cout << "\n The Current Input Image: " << current_image_index + 1 << "/" << ImageNum << " has been read!\n";
			}

			// save parameters into the ".yml" file for image reconstruction
			// 2-kind Parameters of current image include:
			//     * image size: image_h * image_w;
			//     * possible row and column indices of maximun posterior element of each xPatch in currently read image;

				if (!fs_in.isOpened()){
					cerr << "failed to open " << fileStorageName << endl;
					fileStorageHelp();
				}
				else{ // no random shifting, the rox_Idx and col_Idx do not need to be saved
					fs_in << "{:"
						<< "maxPost_row_col_Idx" << v_maxPost_row_column_Idx[current_image_index]
						<< "}";
					if (verbose){
						cout << "v_maxPost_row_column_Idx[" << current_image_index << "] size = " << v_maxPost_row_column_Idx[current_image_index].size() << endl;
						cout << "v_rIdx[" << current_image_index << "] size = " << v_rIdx.size() << endl;
						cout << "v_cIdx[" << current_image_index << "] size = " << v_cIdx.size() << endl;
					}
				}
			

			// delete possible indices of each xPatch in currently read image
			std::vector<int>().swap(v_rIdx);
			std::vector<int>().swap(v_cIdx);

		} // end of each current input image


	// keep the time
	std::cout << "\nThe learning of max-posterior row and column indices finishes now!\n";
	t = clock() - t;
	printf("EM Iteration took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);

	// delete variables and release memory
	// src
	vector<Epitome_Doulbe>().swap(v_eMeanOverVar);
	vector<Epitome_Doulbe>().swap(v_InvVar);
	vector<Epitome_Doulbe>().swap(v_LogVar);
	vector<Epitome_Doulbe>().swap(v_eMean2OverVar);
	//kernel
	vector<Epitome_Doulbe>().swap(v_xPatch);
	vector<Epitome_Doulbe>().swap(v_xxPatch);
	vector<Epitome_Doulbe>().swap(v_patchSize_Ones);
	// Convolution or Correlation Result
	vector<Epitome_Doulbe>().swap(v_eLogVarSum);
	vector<Epitome_Doulbe>().swap(v_eeSum);
	vector<Epitome_Doulbe>().swap(v_double_unnormalized_post);
	vector<Epitome_Doulbe>().swap(v_double_xxSum);
	vector<Epitome_Doulbe>().swap(v_double_xeSum);

	// clear ws
	// get the initial value
	ws.out_src = save_out_src;
	ws.out_kernel = save_out_kernel;
	c_fft::clear_workspace(ws);


	// write the end of the FIleStorage file
	if (!fs_in.isOpened()){
		cerr << "failed to open " << fileStorageName << endl;
	}
	else 
		fs_in << "]";
	fs_in.release();

	// to be continued .......
}