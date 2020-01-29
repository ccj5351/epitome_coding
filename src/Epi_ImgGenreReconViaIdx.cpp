#include "ImageEpitome.hpp"

/*
const std::string DatabaseDir;
const std::string EpitomeResultDir;
const std::string ReconsCompresDir;
const std::string whatkindImgs;
*/
void ImgEpitome::ImgGenreReconViaIdx(
	const vector<string> & v_imageFileNameList // the names of reconstructed images.
	){

	// make a directory for the epitome reconstructed images
	string reconPath = EpitomeResultDir + "/" + nameDifference + to_string(static_cast<long long>(patchSpacing)) 
		+ "_" + to_string(static_cast<long long>(Read_eMean_via_File_Flag)) + "-" + whatkindImgs;
	MakeDir(reconPath);
	int numImages = v_imageFileNameList.size();

	// for each image
	for (int idxImage_itr = 0; idxImage_itr < numImages; ++idxImage_itr){
		// current image size
		int image_h = v_maxPost_row_column_Idx[idxImage_itr][0];
		int image_w = v_maxPost_row_column_Idx[idxImage_itr][1];

		if (verbose){
			cout << "\nThe Current Input Image: " << idxImage_itr + 1 << "/" << numImages << " is being reconstructed.\n"
				<< " >> Its height is " << image_h << " and its width is " << image_w << endl;
		}


		// important parameters for reconstructing an image based on the learned Epitome model
		// before the starting of extracted training patches loop, 
		// please make sure zero-initialized doubles to mat_double_sumQ(..), mat_double_sumQX(..).
		// i.e., clear the matrices used for collecting sufficient statistics.
		
		// zero-initialized doubles
        // must do zero initialization, since due to the following "+=" operation;
		// if no initialized values, the "+=" result can no be controlled.
		cv::Mat mat_double_sumQ_R = cv::Mat::zeros(image_h, image_w, CV_64FC1), 
			mat_double_sumQX_R = cv::Mat::zeros(image_h, image_w, CV_64FC1),
			mat_double_sumQ_G = cv::Mat::zeros(image_h, image_w, CV_64FC1), 
			mat_double_sumQX_G = cv::Mat::zeros(image_h, image_w, CV_64FC1),
			mat_double_sumQ_B = cv::Mat::zeros(image_h, image_w, CV_64FC1), 
			mat_double_sumQX_B = cv::Mat::zeros(image_h, image_w, CV_64FC1),
			mat_double_sumQ_Gray = cv::Mat::zeros(image_h, image_w, CV_64FC1), 
			mat_double_sumQX_Gray = cv::Mat::zeros(image_h, image_w, CV_64FC1);

		///////////////////////////////////////////////////////////////////////////////
		/////////// for each possible xPatchIdx, to generate the most possible patch 
		//////////  based on the learned Epitome, i.e., doing E-Step //////////////////
		///////////////////////////////////////////////////////////////////////////////
		// (row, column) position of the current xPatch

		tuple<vector<int>, vector<int>> t_v_row_col_Idx = getRowColIdxOfOneImg(image_h, image_w);
		vector<int> v_row_Idx = get<0>(t_v_row_col_Idx);
		vector<int> v_col_Idx = get<1>(t_v_row_col_Idx);

		int numRowColIdx = v_maxPost_row_column_Idx[idxImage_itr].size() - 2; // due to the first two elements are image_h and image_w;
		int numXPatches = numRowColIdx;
		int r_size = v_row_Idx.size();
		int c_size = v_col_Idx.size();
		int trainingCounter = 0;
		for (int i = 0, tempIdx = 2 ; i < r_size ; ++i) { // row incides
			for (int j = 0; j < c_size; ++j, ++tempIdx){ // column indices
			int temp_r = v_row_Idx[i];
			int temp_c = v_col_Idx[j];
			int r_maxPosterior = v_maxPost_row_column_Idx[idxImage_itr][tempIdx] / eWidth;
			int c_maxPosterior = v_maxPost_row_column_Idx[idxImage_itr][tempIdx] % eWidth;

			if (verbose & (trainingCounter % (numXPatches / 10) == 0)){
				cout << "  For image " << idxImage_itr + 1 << "/ " << numImages << ":\n";
				std::cout << "    The position of the current xPatch is: " << "row = " << temp_r << "  and column = " << temp_c << endl;
				std::cout << " ---- " << (double(100 * trainingCounter)) / numXPatches << "% Complete\n";
			}

			
			// collect sufficient statistics, taking into consideration that the Epitome wraps around
			// to get the reconstructed Patch from the learned eMean.
			if (isgrayScale){
				for (int i = temp_r; i < temp_r + patchSideLengh; i++) {
					for (int j = temp_c; j < temp_c + patchSideLengh; j++){
						mat_double_sumQ_Gray.at<double>(i, j) += 1;
						// taking into consideration that the Epitome wraps around
						int temp = ((r_maxPosterior + i - temp_r) % eHeight)*eWidth + ((c_maxPosterior + j - temp_c) % eWidth);
						mat_double_sumQX_Gray.at<double>(i, j) += v_eMean_Gray[temp];
					}
				}
			}
			else {
			for (int i = temp_r; i < temp_r + patchSideLengh; i++) {
				for (int j = temp_c; j < temp_c + patchSideLengh; j++){
					// once again:
					// must do zero initialization before this "+=" operation;
					// if no initialized values, the "+=" result can no be controlled.
					mat_double_sumQ_R.at<double>(i, j) += 1;
					mat_double_sumQ_G.at<double>(i, j) += 1;
					mat_double_sumQ_B.at<double>(i, j) += 1;
					// taking into consideration that the Epitome wraps around
					// must do zero initialization before this "+=" operation;
					// if no initialized values, the "+=" result can no be controlled.
					int temp = ((r_maxPosterior + i - temp_r) % eHeight)*eWidth + ((c_maxPosterior + j - temp_c) % eWidth);
					mat_double_sumQX_R.at<double>(i, j) += v_eMean_Red[temp];
					mat_double_sumQX_G.at<double>(i, j) += v_eMean_Gre[temp];
					mat_double_sumQX_B.at<double>(i, j) += v_eMean_Blu[temp];
				}
			}
			}

			trainingCounter++;
		} // end of column
	} // end of row


	// avoid numerical problems
	// to find the nearly-being-zero elements in sumQ
		if (isgrayScale){ /*gray image*/
			for (int r = 0; r < image_h; r++){
				for (int c = 0; c < image_w; c++){
					if (bool numricProblem = mat_double_sumQ_Gray.at<double>(r, c) <= tolerance){
						mat_double_sumQ_Gray.at<double>(r, c) = 1;
						if (verbose){
							cout << "Avoiding numerical problems happens!" << endl;
						}
						// choose a value from eMean
						mat_double_sumQX_Gray.at<double>(r, c) = v_eMean_Gray[(r%eHeight)*eWidth + (c % eWidth)];
					}
				}
			}
		} /*gray image*/
		
		else{ /*color image*/
			for (int r = 0; r < image_h; r++){
				for (int c = 0; c < image_w; c++){
					// Red
					if (bool numricProblem = (mat_double_sumQ_R.at<double>(r, c) <= tolerance)){
						mat_double_sumQ_R.at<double>(r, c) = 1;
						if (verbose){
							cout << "Avoiding numerical problems happens!" << endl;
						}
						// choose a value from eMean
						mat_double_sumQX_R.at<double>(r, c) = v_eMean_Red[(r%eHeight)*eWidth + (c % eWidth)];
					}
					// Green
					if (bool numricProblem = mat_double_sumQ_G.at<double>(r, c) <= tolerance){
						mat_double_sumQ_G.at<double>(r, c) = 1;
						if (verbose){
							cout << "Avoiding numerical problems happens!" << endl;
						}
						// choose a value from eMean
						mat_double_sumQX_G.at<double>(r, c) = v_eMean_Gre[(r%eHeight)*eWidth + (c % eWidth)];
					}
					//Blue
					if (bool numricProblem = mat_double_sumQ_B.at<double>(r, c) <= tolerance){
						mat_double_sumQ_B.at<double>(r, c) = 1;
						if (verbose){
							cout << "Avoiding numerical problems happens!" << endl;
						}
						// choose a value from eMean
						mat_double_sumQX_B.at<double>(r, c) = v_eMean_Blu[(r%eHeight)*eWidth + (c % eWidth)];
					}
				}
			}
		} /*color image*/
	

	// compute the reconstruction
		if (isgrayScale){ /*gray image*/
			Mat temp_ReconImg = Mat(image_h, image_w, CV_64FC1);
			for (int r = 0; r < image_h; r++){
				for (int c = 0; c < image_w; c++){
					temp_ReconImg.at<double>(r, c) = 
						255 * mat_double_sumQX_Gray.at<double>(r, c) / mat_double_sumQ_Gray.at<double>(r, c);
				}
			}

			string filename = reconPath + "/" + nameDifference + v_imageFileNameList[idxImage_itr];
			setFileExtension(filename);
			cv::imwrite(filename, temp_ReconImg);
		} /*gray image*/

		else { /*color image*/
			Mat temp_ReconImg = Mat(image_h, image_w, CV_64FC3);
			for (int r = 0; r < image_h; r++){
				for (int c = 0; c < image_w; c++){
				temp_ReconImg.at<cv::Vec3d>(r, c)[2] = 255 * mat_double_sumQX_R.at<double>(r, c) / mat_double_sumQ_R.at<double>(r, c);
				temp_ReconImg.at<cv::Vec3d>(r, c)[1] = 255 * mat_double_sumQX_G.at<double>(r, c) / mat_double_sumQ_G.at<double>(r, c);
				temp_ReconImg.at<cv::Vec3d>(r, c)[0] = 255 * mat_double_sumQX_B.at<double>(r, c) / mat_double_sumQ_B.at<double>(r, c);
				}
			}
			string filename = reconPath + "/" + nameDifference + v_imageFileNameList[idxImage_itr];
			setFileExtension(filename);
			cv::imwrite(filename, temp_ReconImg);
		}  /*color image*/

	vector<int>().swap(v_row_Idx);
	vector<int>().swap(v_col_Idx);
} // end of each image
}