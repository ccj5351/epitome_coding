#include "ImageEpitome.hpp"

// Note: 
// this function is defined parallel computaiton;
// this means one thread can execute this function
// to output and/or save its result, without conflict
// with other threads.

// This function is designed to:
// get the maximal element p_max_unnormalized_post, 
// and save its corresponding row-column indices into 
// the class member variable -- v_maxPost_row_column_Idx;
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
	// // so this patchIdx should be: 0 <= patchIdx < PatchNum;
	const int  patchIdx,
	const bool IsSaveRowColIdx,
	const int  temp_location, //  temp_location = temp_r * image_w + temp_c;
	const int  image_w,
	const int  current_image_index, // used for saving maxPost_row_column_Idx;
	const cv::Mat & mat_currentInputImage,
	const int  PatchNum, // used for disply the remaining time;
	Epitome_Doulbe * const p_max_unnormalized_post
	){

	if (!mat_currentInputImage.empty()){ /*no empty image*/
		// temp_r and temp_c are row or column indices in the input image;
		int	temp_r = temp_location / image_w;
		int	temp_c = temp_location % image_w;
		int channe = mat_currentInputImage.channels();
		// the pixel number in epitome
		int eNum = eHeight * eWidth;
		// the pixel number in image patch
		int pNum = patchSideLengh * patchSideLengh;

		if (channe == 1){
			/*gray image*/
			vector<Epitome_Doulbe>  v_xPatch_Gray(pNum), v_xxPatch_Gray(pNum);

			vector<Epitome_Doulbe>
				v_xxSum_Gray(eNum), v_xeSum_Gray(eNum),
				v_unnormalized_post(eNum),
				/*Some variables is the summation of all the channels*/
				v_normalized_post(eNum);


			for (int i = temp_r, temp_counter = 0; i < temp_r + patchSideLengh; ++i)
				for (int j = temp_c; j < temp_c + patchSideLengh; ++j, ++temp_counter)
					v_xPatch_Gray[temp_counter] = (double)(mat_currentInputImage.at<uchar>(i, j)) / 255.0f;
			// to get the extracted Patch v_xxPatch.
			for (int temp_counter = 0; temp_counter < pNum; ++temp_counter)
				v_xxPatch_Gray[temp_counter] = v_xPatch_Gray[temp_counter] * v_xPatch_Gray[temp_counter];

			
			///////////////////////////////////////////////////////////////////////////////
			/////////// for xxSum  //////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////
			// correlation
			ws.mode = c_fft::LINEAR_CORRELATION_VALIDSAME;
			update_workspace_few(ws);
			// Gray channel
			c_fft::convolve(ws, &v_InvVar_Gray[0], &v_xxPatch_Gray[0]);
			memcpy(&v_xxSum_Gray[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));

			///////////////////////////////////////////////////////////////////////////////
			/////////// for xeSum  //////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////
			// Gray channel
			c_fft::convolve(ws, &v_eMeanOverVar_Gray[0], &v_xPatch_Gray[0]);
			memcpy(&v_xeSum_Gray[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			//for (int i = 0; i != eNum; ++i)
			//	v_xeSum_Gray[i] = ws.dst_convolution[i];

			// to calculate the unnormalized posterior or likelihood matrix,
			// and find the maximum element wherein and its corresponding index
			* (p_max_unnormalized_post) = DOUBLE_MIN; // to get a minimum value in double

			for (int i = 0; i != eNum; ++i){
				v_unnormalized_post[i] = (-0.5)* (v_eeSum_Gray[i] - 2 * v_xeSum_Gray[i]
					+ v_xxSum_Gray[i] + v_eLogVarSum_Gray[i]);
				if ((*p_max_unnormalized_post) < v_unnormalized_post[i]){
					(*p_max_unnormalized_post) = v_unnormalized_post[i];
					if (IsSaveRowColIdx)
						// +2 , due to the 0-th and 1-th element is the image_h and image_w;
						v_maxPost_row_column_Idx[current_image_index][patchIdx + 2] = i;
				}
			}



			// release memory
			vector<Epitome_Doulbe>().swap(v_xPatch_Gray);
			vector<Epitome_Doulbe>().swap(v_xxPatch_Gray);
			vector<Epitome_Doulbe>().swap(v_xxSum_Gray);
			vector<Epitome_Doulbe>().swap(v_xeSum_Gray);
			vector<Epitome_Doulbe>().swap(v_unnormalized_post);
			return true;

		} /*end of gray image*/
		
		else {/*color image*/
		
			vector<Epitome_Doulbe>  v_xPatch_R(pNum), v_xxPatch_R(pNum),
				v_xPatch_G(pNum), v_xxPatch_G(pNum),
				v_xPatch_B(pNum), v_xxPatch_B(pNum);

			vector<Epitome_Doulbe>
				v_xxSum_R(eNum), v_xeSum_R(eNum),
				v_xxSum_G(eNum), v_xeSum_G(eNum),
				v_xxSum_B(eNum), v_xeSum_B(eNum),
				v_unnormalized_post_R(eNum),
				v_unnormalized_post_G(eNum),
				v_unnormalized_post_B(eNum),
				v_unnormalized_post(eNum),
				/*Some variables is the summation of all the channels*/
				v_normalized_post(eNum);


			// pay attention to Red is channel 2, Green is channel 1, blue is channel 0;
			// Red, channel two;
			for (int i = temp_r, temp_counter = 0; i < temp_r + patchSideLengh; ++i)
				for (int j = temp_c; j < temp_c + patchSideLengh; ++j, ++temp_counter){
						v_xPatch_R[temp_counter] = (double)(mat_currentInputImage.at<cv::Vec3b>(i, j)[2]) / 255.0f;
				}

			// to get the extracted Patch v_xxPatch.
			for (int temp_counter = 0; temp_counter < pNum; ++temp_counter)
				v_xxPatch_R[temp_counter] = v_xPatch_R[temp_counter] * v_xPatch_R[temp_counter];
			// Green, channel one;
			for (int i = temp_r, temp_counter = 0; i < temp_r + patchSideLengh; ++i)
				for (int j = temp_c; j < temp_c + patchSideLengh; ++j, ++temp_counter)
					v_xPatch_G[temp_counter] = (double)(mat_currentInputImage.at<cv::Vec3b>(i, j)[1]) / 255.0f;
			// to get the extracted Patch v_xxPatch.
			for (int temp_counter = 0; temp_counter < pNum; ++temp_counter)
				v_xxPatch_G[temp_counter] = v_xPatch_G[temp_counter] * v_xPatch_G[temp_counter];
			// Blue, channel zero;
			for (int i = temp_r, temp_counter = 0; i < temp_r + patchSideLengh; ++i)
				for (int j = temp_c; j < temp_c + patchSideLengh; ++j, ++temp_counter)
					v_xPatch_B[temp_counter] = (double)(mat_currentInputImage.at<cv::Vec3b>(i, j)[0]) / 255.0f;
			// to get the extracted Patch v_xxPatch.
			for (int temp_counter = 0; temp_counter < pNum; ++temp_counter)
				v_xxPatch_B[temp_counter] = v_xPatch_B[temp_counter] * v_xPatch_B[temp_counter];


			///////////////////////////////////////////////////////////////////////////////
			/////////// for xxSum  //////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////
			// correlation
			ws.mode = c_fft::LINEAR_CORRELATION_VALIDSAME;
			update_workspace_few(ws);
			// Red channel
			c_fft::convolve(ws, &v_InvVar_R[0], &v_xxPatch_R[0]);
			memcpy(&v_xxSum_R[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			//for (int i = 0; i != eNum; ++i)
			//	v_xxSum_R[i] = ws.dst_convolution[i];
			// Green
			c_fft::convolve(ws, &v_InvVar_G[0], &v_xxPatch_G[0]);
			memcpy(&v_xxSum_G[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			//for (int i = 0; i != eNum; ++i)
			//	v_xxSum_G[i] = ws.dst_convolution[i];
			// Blue
			c_fft::convolve(ws, &v_InvVar_B[0], &v_xxPatch_B[0]);
			memcpy(&v_xxSum_B[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			//for (int i = 0; i != eNum; ++i)
			//	v_xxSum_B[i] = ws.dst_convolution[i];


			///////////////////////////////////////////////////////////////////////////////
			/////////// for xeSum  //////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////////////
			// Red channel
			c_fft::convolve(ws, &v_eMeanOverVar_R[0], &v_xPatch_R[0]);
			memcpy(&v_xeSum_R[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			//for (int i = 0; i != eNum; ++i)
			//	v_xeSum_R[i] = ws.dst_convolution[i];
			// Green
			c_fft::convolve(ws, &v_eMeanOverVar_G[0], &v_xPatch_G[0]);
			memcpy(&v_xeSum_G[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			//for (int i = 0; i != eNum; ++i)
			//	v_xeSum_G[i] = ws.dst_convolution[i];
			// Blue
			c_fft::convolve(ws, &v_eMeanOverVar_B[0], &v_xPatch_B[0]);
			memcpy(&v_xeSum_B[0], ws.dst_convolution, eNum * sizeof(Epitome_Doulbe));
			//for (int i = 0; i != eNum; ++i)
			//	v_xeSum_B[i] = ws.dst_convolution[i];

			// to calculate the unnormalized posterior or likelihood matrix,
			// and find the maximum element wherein and its corresponding index
			* (p_max_unnormalized_post) = DOUBLE_MIN; // to get a minimum value in double
			// R channel
			for (int i = 0; i != eNum; ++i)
				v_unnormalized_post_R[i] = (-0.5)* (v_eeSum_R[i] - 2 * v_xeSum_R[i]
				+ v_xxSum_R[i] + v_eLogVarSum_R[i]);

			// G channel
			for (int i = 0; i != eNum; ++i)
				v_unnormalized_post_G[i] = (-0.5)* (v_eeSum_G[i] - 2 * v_xeSum_G[i]
				+ v_xxSum_G[i] + v_eLogVarSum_G[i]);


			// B channel
			for (int i = 0; i != eNum; ++i)
				// int i = r *eWidth + c;
				v_unnormalized_post_B[i] = (-0.5)* (v_eeSum_B[i] - 2 * v_xeSum_B[i]
				+ v_xxSum_B[i] + v_eLogVarSum_B[i]);


			//////////////////////////////////////////////////////////////////////////////////
			// save the corresponding row_column indices of the maximum posterior element 
			// for all the possible xPatch extracted out of the current image
			//////////////////////////////////////////////////////////////////////////////////
			// v_maxPost_row_column_Idx[current_image_index][trainingCounter + 2] = r_maxPosterior* eWidth + c_maxPosterior;

			// do summation over 3 channels;
			for (int i = 0; i < eNum; ++i){
				v_unnormalized_post[i] = v_unnormalized_post_B[i] +
					v_unnormalized_post_G[i] + v_unnormalized_post_R[i];
				if ((*p_max_unnormalized_post) < v_unnormalized_post[i]){
					(*p_max_unnormalized_post) = v_unnormalized_post[i];
					if (IsSaveRowColIdx)
						// +2 , due to the 0-th and 1-th element is the image_h and image_w;
						v_maxPost_row_column_Idx[current_image_index][patchIdx + 2] = i;
				}
			}


		// release memory
		vector<Epitome_Doulbe>().swap(v_xPatch_R);
		vector<Epitome_Doulbe>().swap(v_xxPatch_R);
		vector<Epitome_Doulbe>().swap(v_xxSum_R);
		vector<Epitome_Doulbe>().swap(v_xeSum_R);
		vector<Epitome_Doulbe>().swap(v_unnormalized_post_R);

		vector<Epitome_Doulbe>().swap(v_xPatch_G);
		vector<Epitome_Doulbe>().swap(v_xxPatch_G);
		vector<Epitome_Doulbe>().swap(v_xxSum_G);
		vector<Epitome_Doulbe>().swap(v_xeSum_G);
		vector<Epitome_Doulbe>().swap(v_unnormalized_post_G);

		vector<Epitome_Doulbe>().swap(v_xPatch_B);
		vector<Epitome_Doulbe>().swap(v_xxPatch_B);
		vector<Epitome_Doulbe>().swap(v_xxSum_B);
		vector<Epitome_Doulbe>().swap(v_xeSum_B);
		vector<Epitome_Doulbe>().swap(v_unnormalized_post_B);
		vector<Epitome_Doulbe>().swap(v_unnormalized_post);

		return true;
		} /*end of color image*/
		

} /*end of NOT empty image*/

	else {/*empty image*/
		std::cout << "Error! Image is Empty!\n";
		return false;
	}
}

