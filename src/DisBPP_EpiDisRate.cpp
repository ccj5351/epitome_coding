#include "Dis_Bpp.h"

//*********************************************************************************************************
// this function is used to get the size (unit: Bytes) of all the files (like, images) in some a directory
//*********************************************************************************************************
uintmax_t Distortion_Bpp::getFileSizeFromDir(const string & input_img_dir // input image category directory
	){
	uintmax_t sum_input_size = 0;
	vector<string> filelist;
	GetFileList(input_img_dir, &filelist);
	auto filelistSize = filelist.size();
	for (unsigned int i = 0; i < filelistSize; ++i){
		string input_path = input_img_dir + "/" + filelist[i];
		sum_input_size += file_size(input_path);
	}
	vector<string>().swap(filelist);
	return sum_input_size;
}


//*********************************************************************************************************
// to get rate-distortion data for some compression quality, say, 90, 
// and some encoding and/or quantization parameter setting, say, recon-4_0-Edwin-less;
//*********************************************************************************************************
tuple<double, double, double, double> Distortion_Bpp::getEpiCompressRate_Distortion(
	const string & ReconsCompresDir,
	const uintmax_t & sum_input_size, // input image database file size (Units : Bytes)
	const uintmax_t & epitomeDataFileSize, // epitome Data FileSize, Units = Bytes
	const COMPRE_FLAG_2_EPI_RECON & compres_flag_2_epi_recon
	){
	vector<string> input_filelist;

	GetFileList(DatabaseDir, &input_filelist);
	//vector<tuple<double, double>> v_data(fileListSize, tuple<double, double>());
	//	tuple<0> -- psnr;
	//  tuple<1> -- mse;
	vector<tuple<double, double>> v_data(save_img_number, tuple<double, double>(.0, .0));
	
	uintmax_t sum_errorimg_size = 0;
	// get compressed error/residual image's file size (units : Bytes)
	switch (compres_flag_2_epi_recon){ // e.g., = 0, 1, 2, 3

	case UNIFORM_QUANTIZATION:// # 0
	{
		std::cout << "Do case 0.\n";
		string recon_error_Dir = ReconsCompresDir + "/" + "error-bin";

		// losslessly encode bmp images into .7z files
		string s_7z_Compress_Path = "C:/Users/ccj/Desktop/CPlusPlusProjects_CCJ/7z/7z.exe";
		char * char_7z_Compress_Path = new char[s_7z_Compress_Path.length() + 1];
		std::strcpy(char_7z_Compress_Path, s_7z_Compress_Path.c_str());
		char * input_bin_dir = new char[recon_error_Dir.length() + 1];
		std::strcpy(input_bin_dir, recon_error_Dir.c_str());
		string s_7z_output = recon_error_Dir + ".7z";
		char * output_7z_file = new char[s_7z_output.length() + 1];
		std::strcpy(output_7z_file, s_7z_output.c_str());

		if (!is_File_Exist(output_7z_file)){// if the file does not exist, do compression
			char * p_7z_compress_buffer = new char[MAX_CHAR_NUM_OF_FILES_PATH];
			std::sprintf(p_7z_compress_buffer, "%s a %s %s", char_7z_Compress_Path, output_7z_file, input_bin_dir);
			// compress, from .bmp to .7z;
			std::system(p_7z_compress_buffer);
			delete[] p_7z_compress_buffer;
		}
		else{
			cout << "The file has existed, so do not '.7z' compress!\n";
		}
		// get file size
		sum_errorimg_size = bf::file_size(s_7z_output);
		delete[] char_7z_Compress_Path;
		delete[] input_bin_dir;
		delete[] output_7z_file;
	}
		break;

	case UNIFORM_QUANTIZATION_LOSSLESS_J2K: // # 1
	{
		std::cout << "Do case 1.\n";
		string recon_error_Dir = ReconsCompresDir + "/" + "error-JPEG2K/error-J2K";
		sum_errorimg_size = getFileSizeFromDir(recon_error_Dir);
	}
		break;

	case UNIFORM_QUANTIZATION_LOSSY_JPEG: // # 2
	{
		std::cout << "Do case 2.\n";
		string recon_error_Dir = ReconsCompresDir + "/" + "error-JPEG/error-JPG-lossy";
		sum_errorimg_size = getFileSizeFromDir(recon_error_Dir);
	}
		break;

	case UNIFORM_QUANTIZATION_LOSSY_J2K:// # 3
	{
		std::cout << "Do case 3.\n";
		string recon_error_Dir = ReconsCompresDir + "/" + "error-JPEG2K-lossy/error-J2K-lossy";
		sum_errorimg_size = getFileSizeFromDir(recon_error_Dir);
	}
		break;

	default:
		cout << "Wrong parameters input, and there are only 4 categories of residual compression." << endl;

	} /* end of switch-case */

	unsigned long int NumPixel = 0; // for bbp (bits per pixel)
	// calculate the psnr and mse values of SAVE_IMG_NUM images;
	for (int fileIdx = 0; fileIdx != save_img_number; ++fileIdx){
		
		string input_path = DatabaseDir + "/" + input_filelist[fileIdx];
		
		// image name without file extension
		std::size_t pos = input_filelist[fileIdx].find(".");
		string imageNameNoFileExtension = input_filelist[fileIdx].substr(0, pos);

		// to get PSNR of each image
		Mat_<double> input = imread(input_path, img_read_flag);

		// since the reconstructed images might have been saved into ".jpg" images, due to encodeTypeforErrorImgs = ".jpg";
		// make sure the filename is encodeTypeforReconsImgs, say, ".bmp"
		// setFileExtension(input_filelist[fileIdx], encodeTypeforReconsImgs);
		string recon_path = ReconsCompresDir + "/" + nameDifference + imageNameNoFileExtension + encodeTypeforReconsImgs;
		// cout << recon_path << endl;
		Mat_<double> recon = imread(recon_path, img_read_flag);
		getEpiPSNR_MSE(input, recon, get<0>(v_data[fileIdx]), get<1>(v_data[fileIdx]));

		// to get bits per pixel for each image
		// first we want to save the num of pixels in each image
		 NumPixel += input.rows * input.cols* input.channels();
	} /*end of each image*/

	// to get average psnr 
	// to get average bpp 
	int fileListSize = input_filelist.size();
	double average_psnr = .0, // weighted average PSNR (peak signal-to-noise ratio)
		average_mse = .0; // weighted average MSE (mean square error)
	double average_bpp = .0; // weighted average BPP (bits per pixel)
	double compressionRatio = (double(sum_input_size)) / (double(sum_errorimg_size + epitomeDataFileSize)); // original files' size divided by the compressed files' size
	for (int fileIdx = 0; fileIdx != save_img_number; ++fileIdx){
		average_psnr += get<0>(v_data[fileIdx]);
		average_mse += get<1>(v_data[fileIdx]);
	}
	average_psnr /= save_img_number;
	average_mse /= save_img_number;
	average_bpp = ((double)8 * save_img_number * (epitomeDataFileSize + sum_errorimg_size)) / ((double)(fileListSize * NumPixel));
	// release memory
	vector<tuple<double, double>>().swap(v_data);
	vector<string>().swap(input_filelist);
	return tuple<double, double, double, double>(average_psnr, average_mse, average_bpp, compressionRatio);
}


// /* Note:
//  * This function has a bug that it can deal with files of images to some extend,
//  * More images will not be read and written any more in the form of cv::Mat,
//  * specifically, cv::imread and cv::imwrite funciton fail,
//  * maybe it is due to the memory limit allocated by the OpenCV library for Mat
//	* but it is just a guess, which can not be verified.
void Distortion_Bpp::getEpiCompressRate_DistortionfromDir(
	std::map<std::string,uintmax_t> & m_epitomeDataFileSize, // epitome Data FileSize, Units = Bytes
	const uintmax_t & sum_input_size // input image database file size (Units : Bytes)
	){

	string fileStreamName = ReconsCompresDir + "/" + whatKindofImgs + "-psnr-rmse-bpp-cr.txt";
	// to get the file lists of each input image category,
	// and we should manually guarantee that input images and reconstructed images have the almost same filelist, 
	// just with a constant difference, like "recon_",
	// which makes they can be read via the same string variable.

	vector<string> compressQualityCategories; // e.g., cq-005, cq-010, cq-015, ..., cq-095, cq-100
	vector<string> reconsParamCategories; // e.g., Edwin-less-Epi_0_107_8_4_10_UPDATED_FFTW-4_0-UniQunti-0, etc.

	
	GetDirList(ReconsCompresDir, &compressQualityCategories);
	int compressQualityCategorySize = compressQualityCategories.size();
	int reconsParamCategorySize = 0;
	// for saving the tuples of (psnr, mse, bpp, compression ratio)
	vector<vector<tuple<double, double, double, double>>> v_datas(compressQualityCategorySize);
	ofstream fout_psnr(fileStreamName);

	if (!fout_psnr){
		std::cout << "Can Not Open file " << fileStreamName << endl;
	}


	for (int CQIdx = 0; CQIdx != compressQualityCategorySize; ++CQIdx){ // for each compression quality
		cout << "Compression Quality = " << compressQualityCategories[CQIdx] << endl;
		reconsParamCategories.clear();
		GetDirList(ReconsCompresDir + "/" + compressQualityCategories[CQIdx], &reconsParamCategories);
		reconsParamCategorySize = reconsParamCategories.size();
		fout_psnr << whatKindofImgs << "'s " << compressQualityCategories[CQIdx] + " begins here ----\n";
		v_datas[CQIdx] = vector<tuple<double, double, double, double>>(reconsParamCategorySize, tuple<double, double, double, double>(0., 0., 0., 0.));

		for (int paramCategoryIdx = 0; paramCategoryIdx != reconsParamCategorySize; ++paramCategoryIdx){ // for different cases of epitome size and xPatch side_length

			fout_psnr << "  " + compressQualityCategories[CQIdx] + "--" + reconsParamCategories[paramCategoryIdx] + "--" << endl;

			// for one of the reconsParamCategories
			string s_recon_category = ReconsCompresDir + "/" + compressQualityCategories[CQIdx] + "/" + reconsParamCategories[paramCategoryIdx];

			string tempParams = reconsParamCategories[paramCategoryIdx];
#ifdef _DEBUG
			string sub_tempParams = tempParams.substr(tempParams.length() - 1, 1);
#endif
			COMPRE_FLAG_2_EPI_RECON compres_flag_2_epi_recon;
			// the last character of tempParams, = "0", "1", "2", "3";
			if (tempParams.substr(tempParams.length() - 1, 1) == "0"){
				compres_flag_2_epi_recon = UNIFORM_QUANTIZATION;
			}
			else if (tempParams.substr(tempParams.length() - 1, 1) == "1"){
				compres_flag_2_epi_recon = UNIFORM_QUANTIZATION_LOSSLESS_J2K;
			}
			else if (tempParams.substr(tempParams.length() - 1, 1) == "2"){
				compres_flag_2_epi_recon =  UNIFORM_QUANTIZATION_LOSSY_JPEG;
			}
			else if  (tempParams.substr(tempParams.length() - 1, 1) == "3"){
				compres_flag_2_epi_recon = UNIFORM_QUANTIZATION_LOSSY_J2K;
			}
			else {
				cout << "The last character of tempParams should always be '0', '1', '2', '3'. Error occurs!\n";
			}
			// e.g., extracting the characters "4_0" out of "Edwin-less-Epi_0_107_8_4_10_UPDATED_FFTW-4_0-UniQunti-0";
			string key = reconsParamCategories[paramCategoryIdx].substr(tempParams.length() - 14, 3);// key = 4_0;
			// Maps are associative containers that store elements formed by a combination of a key value and a mapped value, following a specific order.
			tuple<double, double, double, double > t_output = 
				getEpiCompressRate_Distortion(s_recon_category, sum_input_size, m_epitomeDataFileSize[key], compres_flag_2_epi_recon);

			fout_psnr << "    ------ average PSNR(Peak Signal Noise Ratio) = " << get<0>(t_output) << endl;
			fout_psnr << "    ------ average MSE = " << get<1>(t_output) << endl;
			fout_psnr << "    ------ average BPP(Bits Per Pixel) = " << get<2>(t_output) << endl;
			fout_psnr << "    ------ Compression Ratio = " << get<3>(t_output) << endl;
			// save the tuple(PSNR, MSE, bpp, CR)
			v_datas[CQIdx][paramCategoryIdx] = t_output;
			fout_psnr << endl;
		} // end of each recon epitome params
		fout_psnr << endl << endl;
	} // end of each compression quality
	fout_psnr.close();


	// the ".txt" file, which can be used to send to GNUPlot for rate-distortion curve plotting
	// the tuple(PSNR, bpp, CR, RMSE)
	for (int paramCategoryIdx = 0; paramCategoryIdx != reconsParamCategorySize; ++paramCategoryIdx){ // for different cases of epitome size and xPatch side_length
		// save the tuple(PSNR, bpp, CR)
		ofstream fout_psnr_txt(ReconsCompresDir + "/" + reconsParamCategories[paramCategoryIdx] + "-PSNR-RMSE-bpp-cr.txt");
		if (!fout_psnr_txt){
			std::cout << "File Not Opened" << endl;
		}
		else 
			fout_psnr_txt << "# PSNR BPP CR MSE\n"; //for notation of GNUPlot
		for (int CQIdx = 0; CQIdx != compressQualityCategorySize; ++CQIdx){ // for each compression quality
			fout_psnr_txt << get<0>(v_datas[CQIdx][paramCategoryIdx]) << " " // psnr
				<< get<2>(v_datas[CQIdx][paramCategoryIdx]) << " " // bpp
				<< get<3>(v_datas[CQIdx][paramCategoryIdx]) << " " // cr
				<< get<1>(v_datas[CQIdx][paramCategoryIdx]) << endl; // mse
		}  // end of each compression quality
		fout_psnr_txt.close();
	} // end of each recon epitome params	

	// release memory
	vector<string>().swap(compressQualityCategories);
	vector<string>().swap(reconsParamCategories);
	vector<vector<tuple<double, double, double, double>>>().swap(v_datas);
}