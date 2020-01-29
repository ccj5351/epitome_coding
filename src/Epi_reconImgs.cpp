#include "ImageEpitome.hpp"

void ImgEpitome::reconImgsAfterLearning(
	const int & NewpatchSpace, 
	const unsigned int & Read_eMean_Flag){

	// input images' names of the current category
	vector<string> filelist;
	GetFileList(DatabaseDir, &filelist);
	std::string s_numIteration = std::to_string(static_cast<long long>(numIteration));
	
	this ->Read_eMean_via_File_Flag = Read_eMean_Flag;

	////////////////////////////////////
	// Reconstruction
	////////////////////////////////////

	//********************************************************************
	// Step - 1, read necessary data from files
	//********************************************************************
	
	//*******************************************************************
	// read Epitome mean and var from the beforehand saved ".yml" file.
	//*******************************************************************
	
	// read v_eVar from ".yml" file;
	string s_eMeanVar; 
	if (isgrayScale)
		s_eMeanVar = EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMeanVar_EM" + s_numIteration + ".yml";
	else 
		s_eMeanVar = EpitomeResultDir + "/" + whatkindImgs + "_Color_eMeanVar_EM" + s_numIteration + ".yml";
	FileStorage fs_eMeanVar_out(s_eMeanVar, FileStorage::READ);
	if (!fs_eMeanVar_out.isOpened()){
		cerr << "failed to open " << s_eMeanVar << endl;
		fileStorageHelp();
	}
	else {
		if (isgrayScale){
			v_eVar_Gray.clear();
		fs_eMeanVar_out["e_Var_Gray"] >> v_eVar_Gray;
		}
		else{
			v_eVar_Gre.clear();
			v_eVar_Red.clear();
			v_eVar_Blu.clear();
			fs_eMeanVar_out["e_Var_Red"] >> v_eVar_Red;
			fs_eMeanVar_out["e_Var_Gre"] >> v_eVar_Gre;
			fs_eMeanVar_out["e_Var_Blu"] >> v_eVar_Blu;
		}
	}
	
	// read v_Mean from diffeerent files;
	// = 0;

	if (Read_eMean_via_File_Flag == 0){ // read v_eMean form ".yml" file;
		if (!fs_eMeanVar_out.isOpened()){
			cerr << "failed to open " << s_eMeanVar << endl;
			fileStorageHelp();
		}
		else{
			if (isgrayScale){
				v_eMean_Gray.clear();
				fs_eMeanVar_out["e_Mean_Gray"] >> v_eMean_Gray;
			}
			else{
				v_eMean_Red.clear();
				v_eMean_Gre.clear();
				v_eMean_Blu.clear();
				fs_eMeanVar_out["e_Mean_Red"] >> v_eMean_Red;
				fs_eMeanVar_out["e_Mean_Gre"] >> v_eMean_Gre;
				fs_eMeanVar_out["e_Mean_Blu"] >> v_eMean_Blu;
			}

			
		}
		fs_eMeanVar_out.release();
	}


	// = 1;
	else if (Read_eMean_via_File_Flag == 1){ // read v_eMean form ".jpg" file;
		Mat mat_eMean;
		if (isgrayScale){
			s_eMeanVar = EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMean_EM" + s_numIteration + ".jpg";
			Mat mat_eMean =  cv::imread(s_eMeanVar, Epi_Image_Read_Flag);
			if (!mat_eMean.data){// Check for invalid input
				cout << "Could not open or find the image for v_eMean loading from files in Epitome Reconstruction process." << std::endl;
			}
			else{
				for (int r = 0; r < eHeight; r++){
					for (int c = 0; c < eWidth; c++){
						int temp = r*eWidth + c;
						v_eMean_Gray[temp] = mat_eMean.at<uchar>(r, c) / 255.0;
					}
				}
			}
		}
		else { /*color*/
			s_eMeanVar = EpitomeResultDir + "/" + whatkindImgs + "_Color_eMean_EM" + s_numIteration + ".jpg";
			Mat mat_eMean =  cv::imread(s_eMeanVar, Epi_Image_Read_Flag);
			if (!mat_eMean.data){// Check for invalid input
				cout << "Could not open or find the image for v_eMean loading from files in Epitome Reconstruction process." << std::endl;
			}
			else{
				for (int r = 0; r < eHeight; r++){
					for (int c = 0; c < eWidth; c++){
						int temp = r*eWidth + c;
						v_eMean_Red[temp] = mat_eMean.at<cv::Vec3b>(r, c)[2] / 255.0;
						v_eMean_Gre[temp] = mat_eMean.at<cv::Vec3b>(r, c)[1] / 255.0;
						v_eMean_Blu[temp] = mat_eMean.at<cv::Vec3b>(r, c)[0] / 255.0;
					}
				}
			}
		}	 
	}

	else{
		std::cout << "Error occurs for Read_eMean_via_File_Flag. Do nothing!\n";
	}

	//*********************************************************************************
	// Step - 2, learn and save the indices into physical files based on max-posterior
	//*********************************************************************************
	// getMaxPostRowColIdx(filelist, patchSpace);
	paral_getMaxPostRowColIdx(filelist, NewpatchSpace);
	//*********************************************************************
	// read row and column indices from the above saved ".yml" file.
	//*********************************************************************
	this->patchSpacing = NewpatchSpace;
	string s_Idx = EpitomeResultDir + "/" + whatkindImgs + "_" 
		+ to_string(static_cast<long long>(patchSpacing)) + "_" 
		+ to_string(static_cast<long long>(Read_eMean_via_File_Flag)) + "_RowCol_Idx.yml";

	v_maxPost_row_column_Idx.clear();
	v_maxPost_row_column_Idx = read_row_col_idx_form_YML(s_Idx);
	/*
	FileStorage fs_Idx_out(s_Idx, FileStorage::READ);
	if (!fs_Idx_out.isOpened())
	{
		cerr << "failed to open " << s_Idx << endl;
		fileStorageHelp();
	}

	FileNode fn_Idx = fs_Idx_out[whatkindImgs];
	FileNodeIterator it = fn_Idx.begin(), it_end = fn_Idx.end();
	auto numInputImage = fn_Idx.size();
	int temp_idx = 0;

	for (; it != it_end; ++it, temp_idx++){
		v_maxPost_row_column_Idx[temp_idx].clear(); // max-posterior indices read from the beforehand saved files (e.g.".yml")
		(*it)["maxPost_row_col_Idx"] >> v_maxPost_row_column_Idx[temp_idx];
	}

	fs_Idx_out.release();
	*/

	//*********************************************************************
	// Step - 3, the reconstructed images.
	//*********************************************************************
	
	// read image names from the beforehand saved ".yml" file.
	string s_ImaName = EpitomeResultDir + "/" + whatkindImgs + "_imgNames.yml";
	vector<string> v_imageFileNameList;
	FileStorage fs_filelist_out(s_ImaName, FileStorage::READ);
	if (!fs_filelist_out.isOpened())
	{
		cerr << "failed to open " << s_ImaName << endl;
		fileStorageHelp();
	}
	fs_filelist_out[whatkindImgs] >> v_imageFileNameList;
	fs_filelist_out.release();

	// the reconstructed images.
	// to make sure  to update the member variable "patchSpacing" using "NewpatchSpace"
	this->patchSpacing = NewpatchSpace;
	// the reconstructed images.
	ImgGenreReconViaIdx(v_imageFileNameList);

	// delete some vector variable
	vector<string>().swap(v_imageFileNameList);
}

// this function is  made to reconstruct images based on the variant row-column indices
// which are derived from the baseline epitome result via quantization and/or down- or up-sampling
// to the baseline row-column-indices;
void ImgEpitome::reconImgsDirect(
	const int & NewpatchSpace, // patch spacing
	const string & baselineEpitomeResultDir, // to read v_eMeanVar from baseline Epitome result directory;
	const unsigned int & Read_eMean_Flag // 0 (YML), 1 (JPG), 2 (PNG), 3(Binary File), 4(BMP)  
	){

	// input images' names of the current category
	vector<string> filelist;
	GetFileList(DatabaseDir, &filelist);
	std::string s_numIteration = std::to_string(static_cast<long long>(numIteration));

	this->Read_eMean_via_File_Flag = Read_eMean_Flag;

	////////////////////////////////////
	// Reconstruction
	////////////////////////////////////

	//********************************************************************
	// Step - 1, read necessary data from files
	//********************************************************************

	//*******************************************************************
	// read Epitome mean and var from the beforehand saved ".yml" file.
	//*******************************************************************

	// read v_eMeanVar from ".yml" file;
	string s_eMeanVar;
	if (isgrayScale)
		s_eMeanVar = baselineEpitomeResultDir + "/" + whatkindImgs + "_Gray_eMeanVar_EM" + s_numIteration + ".yml";
	else
		s_eMeanVar = baselineEpitomeResultDir + "/" + whatkindImgs + "_Color_eMeanVar_EM" + s_numIteration + ".yml";

	FileStorage fs_eMeanVar_out(s_eMeanVar, FileStorage::READ);


	if (!fs_eMeanVar_out.isOpened()){
		cerr << "failed to open " << s_eMeanVar << endl;
		fileStorageHelp();
	}
	else {
		if (isgrayScale){
			v_eVar_Gray.clear();
			fs_eMeanVar_out["e_Var_Gray"] >> v_eVar_Gray;
		}
		else{
			v_eVar_Gre.clear();
			v_eVar_Red.clear();
			v_eVar_Blu.clear();
			fs_eMeanVar_out["e_Var_Red"] >> v_eVar_Red;
			fs_eMeanVar_out["e_Var_Gre"] >> v_eVar_Gre;
			fs_eMeanVar_out["e_Var_Blu"] >> v_eVar_Blu;
		}
	}

	// read v_Mean from diffeerent files;
	// = 0;

	if (Read_eMean_via_File_Flag == 0){ // read v_eMean form ".yml" file;
		if (!fs_eMeanVar_out.isOpened()){
			cerr << "failed to open " << s_eMeanVar << endl;
			fileStorageHelp();
		}
		else{
			if (isgrayScale){
				v_eMean_Gray.clear();
				fs_eMeanVar_out["e_Mean_Gray"] >> v_eMean_Gray;
			}
			else{
				v_eMean_Red.clear();
				v_eMean_Gre.clear();
				v_eMean_Blu.clear();
				fs_eMeanVar_out["e_Mean_Red"] >> v_eMean_Red;
				fs_eMeanVar_out["e_Mean_Gre"] >> v_eMean_Gre;
				fs_eMeanVar_out["e_Mean_Blu"] >> v_eMean_Blu;
			}


		}
		fs_eMeanVar_out.release();
	}


	// = 1;
	else if (Read_eMean_via_File_Flag == 1){ // read v_eMean form ".jpg" file;
		Mat mat_eMean;
		if (isgrayScale){
			s_eMeanVar = EpitomeResultDir + "/" + whatkindImgs + "_Gray_eMean_EM" + s_numIteration + ".jpg";
			Mat mat_eMean = cv::imread(s_eMeanVar, Epi_Image_Read_Flag);
			if (!mat_eMean.data){// Check for invalid input
				cout << "Could not open or find the image for v_eMean loading from files in Epitome Reconstruction process." << std::endl;
			}
			else{
				for (int r = 0; r < eHeight; r++){
					for (int c = 0; c < eWidth; c++){
						int temp = r*eWidth + c;
						v_eMean_Gray[temp] = mat_eMean.at<uchar>(r, c) / 255.0;
					}
				}
			}
		}
		else { /*color*/
			s_eMeanVar = EpitomeResultDir + "/" + whatkindImgs + "_Color_eMean_EM" + s_numIteration + ".jpg";
			Mat mat_eMean = cv::imread(s_eMeanVar, Epi_Image_Read_Flag);
			if (!mat_eMean.data){// Check for invalid input
				cout << "Could not open or find the image for v_eMean loading from files in Epitome Reconstruction process." << std::endl;
			}
			else{
				for (int r = 0; r < eHeight; r++){
					for (int c = 0; c < eWidth; c++){
						int temp = r*eWidth + c;
						v_eMean_Red[temp] = mat_eMean.at<cv::Vec3b>(r, c)[2] / 255.0;
						v_eMean_Gre[temp] = mat_eMean.at<cv::Vec3b>(r, c)[1] / 255.0;
						v_eMean_Blu[temp] = mat_eMean.at<cv::Vec3b>(r, c)[0] / 255.0;
					}
				}
			}
		}
	}

	else{
		std::cout << "Error occurs for Read_eMean_via_File_Flag. Do nothing!\n";	
	}

	//*********************************************************************
	// Step - 2, read existing quantized and/or to the baseline row-column
	// indices from the ".yml" files.
	//*********************************************************************
	// to make sure  to update the member variable "patchSpacing" using "NewpatchSpace"
	this->patchSpacing = NewpatchSpace;
	string s_Idx = EpitomeResultDir + "/" + whatkindImgs + "_" 
		+ to_string(static_cast<long long>(patchSpacing)) + "_" + 
		to_string(static_cast<long long>(Read_eMean_via_File_Flag)) + "_RowCol_Idx.yml";

	v_maxPost_row_column_Idx.clear();
	v_maxPost_row_column_Idx = read_row_col_idx_form_YML(s_Idx);

	/*
	FileStorage fs_Idx_out(s_Idx, FileStorage::READ);
	if (!fs_Idx_out.isOpened())
	{
		cerr << "failed to open " << s_Idx << endl;
		fileStorageHelp();
	}

	FileNode fn_Idx = fs_Idx_out[whatkindImgs];
	FileNodeIterator it = fn_Idx.begin(), it_end = fn_Idx.end();
	auto numInputImage = fn_Idx.size();
	int temp_idx = 0;

	for (; it != it_end; ++it, temp_idx++){
		v_maxPost_row_column_Idx[temp_idx].clear(); // max-posterior indices read from the beforehand saved files (e.g.".yml")
		(*it)["maxPost_row_col_Idx"] >> v_maxPost_row_column_Idx[temp_idx];
	}

	fs_Idx_out.release();
*/

	//*********************************************************************
	// Step - 3, the reconstructed images.
	//*********************************************************************

	// read image names from the beforehand saved ".yml" file.
	string s_ImaName = baselineEpitomeResultDir + "/" + whatkindImgs + "_imgNames.yml";
	vector<string> v_imageFileNameList;
	FileStorage fs_filelist_out(s_ImaName, FileStorage::READ);
	if (!fs_filelist_out.isOpened())
	{
		cerr << "failed to open " << s_ImaName << endl;
		fileStorageHelp();
	}
	fs_filelist_out[whatkindImgs] >> v_imageFileNameList;
	fs_filelist_out.release();

	// the reconstructed images.
	// to make sure  to update the member variable "patchSpacing" using "NewpatchSpace"
	this->patchSpacing = NewpatchSpace;
	ImgGenreReconViaIdx(v_imageFileNameList);

	// delete some vector variable
	vector<string>().swap(v_imageFileNameList);
}