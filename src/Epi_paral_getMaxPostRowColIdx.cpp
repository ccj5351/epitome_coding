#include "ImageEpitome.hpp"
#include "omp.h"

void ImgEpitome::paral_getMaxPostRowColIdx(
	const vector<string> & v_imageFileNameList, // read input images via their file-names.
	const int & patchSpac){
	// print the consuming time
	clock_t t = clock();
	
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
	char c_numIteration[20];
	sprintf_s(c_numIteration, "%02d", numIteration);
	std::string s_numIteration(c_numIteration);

	cout << "We are going to learn " << whatkindImgs << " category's Epitome."
		<< "\nThis category has " << ImageNum << " images in total.\n";

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

		
		// FFT results
		int Num_FFT_Complex_Outputs = ws.h_fftw * (ws.w_fftw / 2 + 1);

		//#endif


		// EM starts here
		int em_iterator = 1;

			/////////////////////////////////////////////////////////////////////
		    // The definition of all the vector data and FFTW plans /////////////
			///////////////////////////////////////////////////////////////////////

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
					std::cout << "the current input image with " << channels << " channels :"
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


				// the equivalent 1-D indices calculated from the 2-D coordiantes;
				// that is, temp_location_Idx = (temp_r, temp_c) =  temp_r * image_w + temp_c;
				// so, temp_r = temp_location_Idx / image_w ;
				// temp_c = temp_location_Idx % image_w;
				vector <int> v_1D_r_c_Idx(PatchNum);
				for (int r_it = 0, temp = 0; r_it < r_size; ++r_it) {
					for (int c_it = 0; c_it < c_size; ++c_it, ++temp){
						//v_1D_r_c_Idx[r_it * c_size + c_it] = v_rIdx[r_it] * image_w + v_cIdx[c_it];
						v_1D_r_c_Idx[temp] = v_rIdx[r_it] * image_w + v_cIdx[c_it];
					}
				}

				int max_1D_r_c_Idx = r_xEndIdx *image_w + c_xEndIdx;

				// due to the limit of memory space, 
				// we separate the total PatchNum patches, into several chunks,
				// and each chunk has up to MAX_PARALLEL_SIZE patches which will further
				// involve the parallel computation.
				int total_paral_times = std::ceil((double)PatchNum / MAX_PARALLEL_SIZE);
				// but we have to add extra elements to "patchNum", so that it can be divided by MAX_PARALLEL_SIZE
				// with no remainder.
				int oldPatchNum = PatchNum; // save in case of later use
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
				cout << "After random value added, max_1D_r_c_Idx = " << *max_element(v_1D_r_c_Idx.begin(), v_1D_r_c_Idx.end()) << endl;
#endif // DEBUG

				int nThreads = 1;

				// the corresponding row_column indices of the maximum posterior element 
				// for all the possible xPatch extracted out of the current image
				
					// due to the last time result left within this variable,
					// we do not need it;
					/*vector.clear() : Clear content
					Removes all elements from the vector (which are destroyed), 
					leaving the container with a size of 0.
					A reallocation is not guaranteed to happen, 
					and the vector capacity is not guaranteed to change due to calling this function. 
					A typical alternative that forces a reallocation is to use swap:
					vector<T>().swap(x);  // clear x reallocating
					*/
				if (v_maxPost_row_column_Idx[current_image_index].size() != 0)
					v_maxPost_row_column_Idx[current_image_index].clear();
					// we will save the new maxPost_row_column_Idx
					v_maxPost_row_column_Idx[current_image_index] = vector<int>(PatchNum + 2, 0);	
				v_maxPost_row_column_Idx[current_image_index][0] = image_h;
				v_maxPost_row_column_Idx[current_image_index][1] = image_w;
				

				// for saving the max_unnormalized_post for each image patch out of one image;
				v_max_unnormalized_post = vector<Epitome_Doulbe>(PatchNum);

				//****************************************************************************
				//****************************************************************************
				//parallel calculation begins here for each xPatch(R, G, B, and/or Gray)
				//****************************************************************************
				//****************************************************************************

				long int trainingCounter = 0; // counting the number of training patches for showing the remaining time;


				bool IsSaveRowColIdx = true;
				// due to the limit of memory space, 
				// we separate the total PatchNum patches, into several chunks,
				// and each chunk has up to MAX_PARALLEL_SIZE patches which will further
				// involve the parallel computation.

				/*several parallel computation*/
				for (int paral_time = 0; paral_time < total_paral_times; ++paral_time){
					// par attention to those calculation,
					// _1D_Idx_iter = idx_start; _1D_Idx_iter < idx_end; ++_1D_Idx_iter
					// they can guarantee that _1D_Idx_iter can be every number in the range [0, PatchNum-1]
					int idx_start = paral_time* MAX_PARALLEL_SIZE;
					int idx_end = (paral_time + 1)* MAX_PARALLEL_SIZE;

					omp_set_num_threads(NUM_THREADS);
					// omp_set_num_threads(1);

#pragma omp parallel default (none) shared (v_ws, idx_start, idx_end, nThreads, trainingCounter, v_1D_r_c_Idx, IsSaveRowColIdx, current_image_index, PatchNum, ImageNum, image_w, mat_currentInputImage)
					{

#pragma omp master
						nThreads = omp_get_num_threads();
#ifdef _DEBUG
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

						for (int _1D_Idx_iter = idx_start; _1D_Idx_iter < idx_end; ++_1D_Idx_iter){
							// paral_time

							bool flag = (_1D_Idx_iter == 998) ? true : false;
							int temp = 0;
							if (get_xPatch_unnormalPosterior_RGB_or_Gray_Channel(
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
								&v_max_unnormalized_post[_1D_Idx_iter]
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
					omp_set_num_threads(1);
					// and check whether it has been successfully set;
					nThreads = omp_get_num_threads();
#ifdef _DEBUG
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
			//*******************************************************
			// now save all the indices into yml files;
			//*******************************************************
			char c_temp[20];
			sprintf_s(c_temp, "%d_%d", patchSpacing, Read_eMean_via_File_Flag);
			std::string s_patchSpacing(c_temp);

			std::string fileStorageName = EpitomeResultDir + "/"
				+ whatkindImgs + "_" + c_temp + "_RowCol_Idx.yml";

			write_row_col_idx_2_YML(fileStorageName, ImageNum);
			/*
			FileStorage fs_in(fileStorageName, FileStorage::WRITE);
			if (!fs_in.isOpened())
			{
				cerr << "failed to open " << fileStorageName << endl;
				fileStorageHelp();
			}
			else // just write the beginning of the FileStorage file.
				fs_in << whatkindImgs << "[:";

			for (int current_image_index = 0; current_image_index != ImageNum; ++current_image_index){
				
				if (!fs_in.isOpened()){
					cerr << "failed to open " << fileStorageName << endl;
					fileStorageHelp();
				}
				else{ // no random shifting, the rox_Idx and col_Idx do not need to be saved
					fs_in << "{:"
						<< "maxPost_row_col_Idx" << v_maxPost_row_column_Idx[current_image_index]
						<< "}";
					if (verbose)
						cout << "v_maxPost_row_column_Idx[" << current_image_index << "] size = " << v_maxPost_row_column_Idx[current_image_index].size() << endl;
					}
			}
		
			// write the end of the FIleStorage file
			if (!fs_in.isOpened()){
				cerr << "failed to open " << fileStorageName << endl;
			}
			else
				fs_in << "]";
			fs_in.release();
			*/

			// keep the time
			std::cout << "\nThe learning of max-posterior row and column indices finishes now!\n";
			t = clock() - t;
			printf("EM Iteration took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);

			// clear ws
			// get the initial value
			ws.out_src = save_out_src;
			ws.out_kernel = save_out_kernel;
			c_fft::clear_workspace(ws);

		// to be continued .......
		// to be continued .......
		// to be continued .......
}