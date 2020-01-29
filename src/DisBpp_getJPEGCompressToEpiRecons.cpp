#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include "Dis_Bpp.h"
#include "combine4BitsImg.h"
#include "quantizer.h"
#include "bitstream.h"


void Distortion_Bpp::getJPEGCompressToEpiRecons(bool & IsErrorImghistgram){

	// to get the file lists of each input image category,
	// and we should manually guarantee that input images and reconstructed images have the almost same filelist, 
	// just with a constant difference, like "recon-",
	// which makes they can be read via the same string variable.
	vector<string> reconsParamCategories, filelist;
	GetDirList(EpitomeResultDir, &reconsParamCategories);



	// for JPEG compression quality
	// For JPEG, different compression quality
	// it can be a quality ( CV_IMWRITE_JPEG_QUALITY ) from 0 to 100 (the higher is the better). Default value is 95.			
	// vector<int> v_jpeg100 = { CV_IMWRITE_JPEG_QUALITY, 100 };

	vector<int> v_jpegParams(2);
	v_jpegParams[0] = CV_IMWRITE_JPEG_QUALITY;
	v_jpegParams[1] = jpeg_quality;
	
	// /* check the string is 1-, 2-, or 3-character length, and to make them be 3-character length, if less than 3.
	//    i.e., change '5' to '005', and '15' to '015'; 
	//    Usage: std::string str ("Test string");
	//    std::cout << "The size of str is " << str.length() << " bytes.\n";
	// */
	int imgReadFlag = 0; // for gray images

	// number char control
	char c_jpeg_quality[10];
	sprintf(c_jpeg_quality, "%03d", jpeg_quality);
	string s_jpeg_q(c_jpeg_quality);
	/*
	if (s_jpeg_q.length() == 1){
		s_jpeg_q = "00" + s_jpeg_q;
	}

	if (s_jpeg_q.length() == 2){
		s_jpeg_q = "0" + s_jpeg_q;
	}
	*/

	//  char * s_jpeg_q; 
	// sprintf(s_jpeg_q, "%03d", jpeg_quality);


	string fileStreamForErrorImgSize = ReconsCompresDir + "/" + "errorImgSize-cq-" + s_jpeg_q + ".txt";
	ofstream fout_errorImgSize(fileStreamForErrorImgSize, std::ofstream::app);
	if (!fout_errorImgSize){
		std::cout << "File Not Opened" << endl;
	}
	s_jpeg_q = "cq-" + s_jpeg_q;

	MakeDir(ReconsCompresDir + "/" + s_jpeg_q);

	
	int reconsParamCategorySize = reconsParamCategories.size();

	for (int paramCategoryIdx = 0; paramCategoryIdx != reconsParamCategorySize; ++paramCategoryIdx){ // for different cases of epitome size and xPatch side_length

		fout_errorImgSize << s_jpeg_q + " -- " + reconsParamCategories[paramCategoryIdx] + " -- " + s_quantizeType << " - error images size :: " << endl;

		// make directories for saving the final compressed reconstruction images
		string tempDir = ReconsCompresDir + "/" + s_jpeg_q + "/" + reconsParamCategories[paramCategoryIdx]
		+ "-" + "errHistThres-" + to_string(static_cast<long long>(errHistThres)) + "-" + s_quantizeType;
		MakeDir(tempDir);
		// make a directory for the error images, whose sizes are needed for the following rate-distortion calculation.
		string nameErrorImgsDir = tempDir + "/" + nameDifference + "error";
		// here we make a separate directory for those recon_error images, just because it is easy to measure their total sizes,
		// when they are in the same directory.
		MakeDir(nameErrorImgsDir);
		// get image names in each category
		GetFileList(DatabaseDir, &filelist);
		int fileListSize = filelist.size();
		for (int fileIdx = 0; fileIdx != fileListSize; ++fileIdx){
			// read input image
			const Mat_<double> input = imread(DatabaseDir + "/" + filelist[fileIdx], imgReadFlag); // read gray images
			// read reconstructed image
			// first check the file name extension
			setFileExtension(filelist[fileIdx], encodeTypeforReconsImgs);
			vector<string> patchSpacingDir;
			GetDirList(EpitomeResultDir + "/" + reconsParamCategories[paramCategoryIdx], &patchSpacingDir);
			const Mat_<double> reconsInput = imread(EpitomeResultDir + "/" + reconsParamCategories[paramCategoryIdx] + "/" + patchSpacingDir[0] +"/"
				                             + nameDifference + filelist[fileIdx], imgReadFlag);

			// do some normalization for the following JPEG compression and saving
			Mat_<double> error = input - reconsInput;
			
			// keep the max and min values of error, before normalization.
			double max_intensity, min_intensity;
			minMaxLoc(error, &min_intensity, &max_intensity);
			
			//********************************************
			// some variables for uniform quantization
			//********************************************
			// keep the unnormalized value for later use
			Mat_<double> unnormalized_error = input - reconsInput;
			// here we use uchar, since we will guarantee that the quantized value level is equal to or less than 256.
			Mat_<uchar> quantized_unnormalized_error = Mat::zeros(input.rows, input.cols, CV_8UC1);
			int nLevels;
			// Function : FILE * fopen ( const char * filename, const char * mode );
			string error_binary_header_file_path;
			const char * char_error_binary_header_file_path;
			FILE * p_write_header_file = NULL;
			string error_binary_file_path;
			const char * char_error_binary_file_path;
			FILE * p_write_file = NULL;
			BitStream bs(p_write_header_file);
			
			// do normalization
			// intensity 0 mapping to 127, i.e., 0 ----> 127 (median)
			// after normalization, error belonging to the range [0, 255]
			double maxAbs = (abs(max_intensity)) > (abs(min_intensity)) ? abs(max_intensity) : abs(min_intensity);
			int maxItensity = pow(2.0, errorImgBit) - 1; // e.g., == 255;
			int midItensity = pow(2.0, errorImgBit - 1) - 1; // e.g. == 127;
			error = (maxItensity / (2 * maxAbs)) * error + midItensity;
			for (int i = 0, rowNum = error.rows; i != rowNum; ++i)
				for (int j = 0, colNum = error.cols; j != colNum; ++j){
				error.at<double>(i, j) = std::max<double>(error.at<double>(i, j), 0.0);// >= 0;
				error.at<double>(i, j) = std::min<double>(error.at<double>(i, j), 255.0);// <= 255
				}

#ifdef _DEBUG
			cout << "Normalized_error :\n";
			display_n_by_n_elements(error, 10);
			// namedWindow("Normalized_error");
			// imshow("Normalized_error", error);

#endif
			// the final compressed recons
			Mat recons;
			// /* to easily specify the compression factor when compressing images on openCV without having to declare a dummy vector.
			// * As suggested, activating C++0x allows me to pass a vector explicitly defined inline to the function. This solved the issue. 
			// * std::vector<int>({1,2}) inline 
			// * If your compiler supports C++11, you can simply do: std::vector<int> v = { 1, 2, 3, 4 };
			// */



			// Firstly, to check the file extension is ".jpg" or not
			// For JPEG, it can be a quality ( CV_IMWRITE_JPEG_QUALITY ) from 0 to 100 (the higher is the better). Default value is 95.
			// make sure the filename is imencodeType, say, ".jpg"
			setFileExtension(filelist[fileIdx], encodeTypeforErrorImgs);
			// std::size_t pos = filelist[fileIdx].find(".");      // e.g., position of ".jpg" in nameTemp
			// filelist[fileIdx] = filelist[fileIdx].substr(0, pos) + encodeTypeforErrorImgs;     // get beginning to ".jpg"

			string error_save_path = nameErrorImgsDir + "/" + nameDifference + "error-" + filelist[fileIdx];
			
			// combine_4_Bits
			bool is_combine4BitsImg = false;
			combine_4_Bit_Img cm4bits;
			cm4bits.initil(); // dynamic build uchar array via "new" operator;

			if (IsErrorImghistgram){
				string pathHist = tempDir + "/" + "errorHistogram";
				MakeDir(pathHist);
				const string foutHistVal = pathHist + "/" + "error-histogram.txt";
#ifdef _DEBUG
				cout << "foutHistVal File Name = " << foutHistVal << "\n";
#endif
				pathHist = pathHist + "/" + "err-hist-" + filelist[fileIdx];

				/// Establish the number of bins
				int histSize = 256;
				int hist_w = 512; //  or 1024, width of histogram drawing
				int hist_h = 400; // height of histogram drawing
				bool  IsCurveTypeHist = false;
				tuple<int, int, int, double, double> max_min_filtered_intensity = 
					getGrayImgHistogram(error, histSize, hist_w, hist_h, IsCurveTypeHist, pathHist, foutHistVal);
				double min_FilteredIntensity = (double)get<0>(max_min_filtered_intensity); // e.g., histogram threshold value == 0.5%
				double med_FilteredIntensity = (double)get<1>(max_min_filtered_intensity);
				double max_FilteredIntensity = (double)get<2>(max_min_filtered_intensity);
				double mean_Intensity = get<3>(max_min_filtered_intensity); // mean
				double sd_Intensity = get<4>(max_min_filtered_intensity); // standard deviation
				
				// 
				double deltaError = maxItensity * errHistThres / 100;
				double errorUpperBound = med_FilteredIntensity + deltaError;
				double errorLoweround = med_FilteredIntensity - deltaError;
 

				Mat_<uchar> combine_error_bit_plane, error_bit_plane, error_unchar;
				switch (quantizeType){
				case SET2MEDIAN: // quantization method 1;
					for (int i = 0, rowNum = error.rows; i != rowNum; ++i){
						for (int j = 0, colNum = error.cols; j != colNum; ++j){
							if (error.at<double>(i, j) >= errorLoweround | error.at<double>(i, j) <= errorUpperBound)
								error.at<double>(i, j) = med_FilteredIntensity;
						}
					}
					imwrite(error_save_path, error, v_jpegParams);
					break;

				case BETWEEN_MIN_MAX:
					for (int i = 0, rowNum = error.rows; i != rowNum; ++i){
						for (int j = 0, colNum = error.cols; j != colNum; ++j){
							error.at<double>(i, j) = std::max<double>(error.at<double>(i, j), min_FilteredIntensity);
							error.at<double>(i, j) = std::min<double>(error.at<double>(i, j), max_FilteredIntensity);
						}
					}
					imwrite(error_save_path, error, v_jpegParams);
					break;
				
				case NOTHING:
					// do nothing to error image
					imwrite(error_save_path, error, v_jpegParams);
					break;

				// The Visual Studio automatically defines _DEBUG symbol for Debug builds (and NDEBUG for non-debug builds).
				case BIT_PLANE_SLICING_JPEG:
					// do compression to error via bit-plane-slicing
					error.convertTo(error_unchar, CV_8UC1);
					// error_bit_plane = recon_via_2_MSB_BitPlane(error_unchar);
					error_bit_plane = recon_via_4_MSB_BitPlane(error_unchar);
                    #ifdef _DEBUG
					namedWindow("error_bit_plane");
					imshow("error_bit_plane", error_bit_plane);
					cout << "reconstruction of error via bit_plane_slicing :\n";
					display_n_by_n_uchar_elements(error_bit_plane, 10);
					// setFileExtension(filelist[fileIdx], ".bmp");
					// nameTemp2 = nameErrorImgsDir + "/" + nameDifference + "error-" + filelist[fileIdx];
					// cout << nameTemp2 << endl;
                    #endif
					imwrite(error_save_path, error_bit_plane, v_jpegParams);

					break;
				
				case BIT_PLANE_SLICING_COMBINE4BITSIMG_JPEG:
					// do compression to error via bit-plane-slicing
					error.convertTo(error_unchar, CV_8UC1);
					error_bit_plane = recon_via_4_MSB_BitPlane(error_unchar); 
					combine_error_bit_plane = cm4bits.Combine(error_bit_plane);
					is_combine4BitsImg = true;					
                    #ifdef _DEBUG
					cout << "combined error via bit_plane_slicing :\n";
					display_n_by_n_uchar_elements(combine_error_bit_plane, 10);
					namedWindow("combine_error_bit_plane");
					imshow("combine_error_bit_plane", combine_error_bit_plane);
                    #endif
					// do JPEG compression
					imwrite(error_save_path, combine_error_bit_plane, v_jpegParams);
					break;

				case UNIFORM_QUANTIZATION_JPEG:
				{
				
					// letting nLevels = 256 actually takes the function of normalization of error to [0, 255] range.
					nLevels = 256; // i.e., 0 - 255
					UniformQuant uniquant_JPEG(max_intensity, min_intensity, nLevels);
					UniformQuantize(unnormalized_error, quantized_unnormalized_error, uniquant_JPEG);

					// do JPEG compression
					imwrite(error_save_path, quantized_unnormalized_error, v_jpegParams);
					break;
				}

				case UNIFORM_QUANTIZATION_BINARY_FILE:
				{
					//****************************************************
					// do uniform quantization here!
					//****************************************************
					
					if (jpeg_quality == 0){
						nLevels = int(max_intensity - min_intensity);
					}
					else
						nLevels = (int)(jpeg_quality * (max_intensity - min_intensity)/100); 
					nLevels = nLevels > 256 ? 256 : nLevels;
					// double log2 (double x);
					// double log2 (T x);  // additional overloads for integral types
					// #include <math.h>       /* ceil */
					// double ceil (double x);
					// ceil () -- returns the smallest integral value that is not less than x (as a floating-point value).
					int bit_for_quantized_error = (int) ceil(log2(nLevels));
					UniformQuant uniQuant_BINARY_FILE(max_intensity, min_intensity, nLevels);
					// Function : FILE * fopen ( const char * filename, const char * mode );
					setFileExtension(error_save_path, ".bin");
					//**************************
					// for header binary file
					//**************************
					error_binary_header_file_path = nameErrorImgsDir + "/" + nameDifference + "error-header-" + filelist[fileIdx];
					char_error_binary_header_file_path = error_binary_header_file_path.c_str();
					p_write_header_file = fopen(char_error_binary_header_file_path, "wb");
				

					if (!p_write_header_file){ // if NULL
						cout << "Cannot open the file of " << char_error_binary_header_file_path << endl;
					}

					else{
						bs.SetFilePointer(p_write_header_file);
						uniQuant_BINARY_FILE.WriteHeader(bs);
					}

					// for body binary file
					error_binary_file_path = error_save_path;
					char_error_binary_file_path = error_binary_file_path.c_str();
					p_write_file = fopen(char_error_binary_file_path, "wb");
					if (!p_write_file){ // if NULL
						cout << "Cannot open the file of " << char_error_binary_file_path << endl;
					}
					else { // write bit-stream
						bs.SetFilePointer(p_write_file);
						UniformQuantize(unnormalized_error, quantized_unnormalized_error, uniQuant_BINARY_FILE);
						for (int i = 0, i_size = quantized_unnormalized_error.rows; i != i_size; ++i){
							for (int j = 0, j_size = quantized_unnormalized_error.cols; j != j_size; ++j){
								bs.WriteBits(quantized_unnormalized_error.at<uchar>(i, j), bit_for_quantized_error); // bit_for_quantized_error
							}
						}
#ifdef _DEBUG
					cout << "before bs.Flush(), numBitsWritten =" << bs.numBitsWritten;
#endif
						
						bs.Flush();
						cout << "After bs.Flush(), numBitsWritten =" << bs.numBitsWritten;

					}
					
					break;
				}
				default:
					cout << "No quantization has been applied to the residual image.\n";
					break;
				}
 
				Mat_<double> errorTemp;
				
				if (is_combine4BitsImg){
					Mat combine_temp = imread(error_save_path, imgReadFlag);
					Mat errorTemp2 = cm4bits.UnCombine(combine_temp);
					// rowRange ¨C Range of the m rows to take. 
					// As usual, the range start is inclusive and the range end is exclusive. 
					// Use Range::all() to take all the rows.
					// colRange ¨C Range of the m columns to take.Use Range::all() to take all the columns.
					errorTemp = Mat::Mat(errorTemp2, Range::all(), Range::Range(0, error.cols));
#ifdef _DEBUG
					cout << "de-combined error via bit_plane_slicing :\n";
					display_n_by_n_elements(errorTemp, 10);
					namedWindow("errorTempDeCombinedBits");
					imshow("errorTempDeCombinedBits", errorTemp);
					
#endif
				}
				else{
					if (quantizeType == UNIFORM_QUANTIZATION_BINARY_FILE){

					}
					else if (quantizeType == UNIFORM_QUANTIZATION_JPEG){

					}
					Mat_<uchar> errorTemp_2 = imread(error_save_path, imgReadFlag);
					errorTemp = Mat_<double>::zeros(errorTemp_2.size());
					for (int i = 0, size_i = errorTemp_2.rows; i != size_i; ++ i){
						for (int j = 0, size_j = errorTemp_2.cols; j != size_j; ++j){
							errorTemp.at<double>(i, j) = (double)errorTemp_2.at<uchar>(i, j);
						}
					}
#ifdef _DEBUG
					cout << "uchar error read from files:\n";
					display_n_by_n_uchar_elements(errorTemp_2, 10);
					cout << "doubled error :\n";
					display_n_by_n_elements(errorTemp, 10);
					cout << errorTemp.at<double>(9, 9) << endl;
					namedWindow("errorTempNoCombinedBits");
					imshow("errorTempNoCombinedBits", errorTemp);					
#endif
				}
				
				cm4bits.clear();
				//		imencode(imencodeType, error, buf, v_jpegParams);
				// read or decode the just saved normalized error image
				// Mat_<double> errorTemp = imdecode(buf, imgReadFlag); // read gray-level images

				// 

				// De-normalize the normalized error
				errorTemp = (errorTemp - midItensity) * 2 * maxAbs / maxItensity;
				recons = Mat(reconsInput + errorTemp);

				// make sure the filename is "*.bmp", since we want to get a lossless recons for visualization
				setFileExtension(filelist[fileIdx], encodeTypeforReconsImgs);
				// imread() for PNG images, it can be the compression level ( CV_IMWRITE_PNG_COMPRESSION ) from 0 to 9. 
				// A higher value means a smaller size and longer compression time. Default value is 3.
				// vector<int> v_png_0 = { CV_IMWRITE_PNG_COMPRESSION, 0 };
				string nameTemp = tempDir + "/" + nameDifference + filelist[fileIdx];
				imwrite(nameTemp, recons);

			} //end of histogram

			//std::vector::clear
			// void clear(); // Clear content, i.e., Removes all elements from the vector (which are destroyed), leaving the container with a size of 0.
			// Note: A reallocation is not guaranteed to happen, and the vector capacity is not guaranteed to change due to calling this function.
			// A typical alternative that forces a reallocation is to use swap:
			// vector<T>().swap(x);   // clear x reallocating
			filelist.clear(); // clear the filelist variable, otherwise, its size will increase, since the new elements, created by the next for-loop, 
			//will not override the old elements, but be added following the old ones.
			//cout << endl;

			// get the file-size of error images, (unit: Bytes) 
			uintmax_t tempSizeBytes = getFileSizeFromDir(nameErrorImgsDir); // input image category directory
			fout_errorImgSize << whatKindofImgs << " = " << tempSizeBytes << "  Bytes\n";

			fout_errorImgSize << endl << endl;
		} //end of "filelist" for-loop
		fout_errorImgSize.close();
	} // end of "reconsParamCategorySize" for-loop
}