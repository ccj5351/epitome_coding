#define _epitome_fftw_convolve_time_test_
#ifndef  _epitome_fftw_convolve_time_test_
#include "allHeaders.hpp"

int main(int argc, char * argv[]){
	clock_t t;
	t = clock();
	std::cout << " Calculating...\n";

	// ***********************************
	// ******* general parameters ********
	// ***********************************
	// for filepath, please use "/" or "\\", not "\", e.g.  "E:\ComputerVisionCCJ\Epitome\imageDatabase\train", 
	// or  "E:\\ComputerVisionCCJ\\Epitome\\imageDatabase";
	string DatabaseDir = "C:/Users/ccj/imageDatabase/bmp-stevens-sub-images/Edwin-less"; //argv[1];
	string EpitomeResultDir = "C:/Users/ccj/imageDatabase/bmp-images-epi-result/epi-6_Try_epi-sub-Edwin-less"; // argv[2];
	string ReconsCompresDir = "C:/Users/ccj/imageDatabase/bmp-images-compress-result/epi-6_Try_compre-sub-Edwin-less-d";//argv[3];
	string GetJPEGCompressDir = ReconsCompresDir;//argv[3];
	bool display= false; // argv[4]
	string whatkindImgs = "Edwin-less"; //argv[5];
	const int newPatchSpacingNum = 5;//  atoi(argv[6]); // e.g., recon-4, recon-8, different values of newPatchSpacing;
	const string s_newPatchSpacing = "4-5-6-7-8"; // argv[7];
	int newPatchSpacing = 0;

	// gain some important parameters from the input parameters of the main function, i.e., argv[i] ;
	int mainFuncParamsIdx = 8; // define some a number for i, just for convenience when some new parameters are added.
	bool randomshifting = false; // make a random shift form the regular patches' positions
	 
	int epitomeWidth = 200;
	int epitomeHeight = epitomeWidth;
	int patchSideLengh = 8;
	int patchSpacing = 4;
	vector<string > inputImgsList;
	GetFileList(DatabaseDir, &inputImgsList);
	int numImgs = inputImgsList.size();
	int numIteration = 10;
	std::string s_numIteration = "10";
	std::string nameDifference = "recon-";
	std::string imgEpi_EncodeType = ".bmp";
	int overhead =0; // file size units = Bytes, e.g., epitome Data FileSize, Units = Bytes
	std::string s_flag_epi_fftw = "UPDATED_FFTW";

	// display some important input parameters of main function
	std::cout << " //*****************************************************" << endl;
	std::cout << " //***display input parameters of main function*********" << endl;
	std::cout << " //*****************************************************" << endl;
	std::cout << " Input Image Dir = " << DatabaseDir << endl; // #1
	std::cout << " Epitome Recons Dir = " << EpitomeResultDir << endl; // #2
	std::cout << " Recons Compression Dir = " << ReconsCompresDir << endl; // #3
	if (display != false)
		std::cout << " display = " << "true" << endl; // #4
	else
		std::cout << " display = " << "false" << endl; // #4
	std::cout << " whatkindImgs = " << whatkindImgs << endl // #5
		<< " newPatchSpacingNum = " << newPatchSpacingNum << endl // #6
		<< " s_newPatchSpacing = " << s_newPatchSpacing << endl; // #7



#ifdef _DEBUG
	std::cout << " recon PatchSpacing s_newPatchSpacing = " << s_newPatchSpacing << endl;
#endif // _DEBUG

	if (randomshifting != false)
		std::cout << " randomshifting = " << "true" << endl; // +0
	else
		std::cout << " randomshifting = " << "false" << endl;// +0

	std::cout << " epitomeWidth =  epitomeHeight = " << epitomeHeight << endl; // +1
	std::cout << " patchSideLengh = " << patchSideLengh << endl; // +2
	std::cout << " patchSpacing = " << patchSpacing << endl; // +3
	std::cout << " EM Iteration  = " << numIteration << endl; // +4
	std::cout << " nameDifference  = " << nameDifference << endl; // +5
	std::cout << " Image Encoding Type  = " << imgEpi_EncodeType << endl; // +6
	std::cout << " Overhead Size = " << overhead << endl; // +7
	std::cout << " UPDATED_FFTW_flag = " << s_flag_epi_fftw << endl; // +8

	std::cout << " //*****************************************************" << endl;
	std::cout << " //*****************************************************" << endl;


#define step_reconstrcuction
#ifndef  step_reconstrcuction
	// ***************************************************
	// ******* initialize GrayImgEpitome Object **********
	// ***************************************************
	// make some necessary directories 
	// epitome reconstruction dir
	MakeDir(EpitomeResultDir);
	EpitomeResultDir = EpitomeResultDir + "/" + whatkindImgs + "-Epi_" + s_randomshifting + "_" + std::to_string(epitomeWidth) + "_"
		+ std::to_string(patchSideLengh) + "_" + std::to_string(patchSpacing) + "_" + std::to_string(numIteration) + "_" + s_flag_epi_fftw;
	MakeDir(EpitomeResultDir);
	ImgEpitome * img1 = new ImgEpitome(epitomeWidth, epitomeHeight, numImgs, patchSideLengh, patchSpacing, numIteration, display,
		randomshifting, DatabaseDir, EpitomeResultDir, ReconsCompresDir, whatkindImgs, nameDifference, imgEpi_EncodeType, overhead);

	string tempImgPath = DatabaseDir + "/" + inputImgsList[0];
	int  imgReadFlag = 0; // 0 = gray ; 1 = color
	Mat mat_currentInputImage = imread(tempImgPath, imgReadFlag);

	bool IsConstantValueInitial = false;
	double constantVal = 0.15;
	img1->initialImgEpitome(mat_currentInputImage, IsConstantValueInitial, constantVal);

	////////////////////////////////////
	// Step - Learn Epitome
	////////////////////////////////////
	EpitomeVariableName compres_flag;
	if (s_flag_epi_fftw == "ORIGINAL_FFTW")
		compres_flag = ORIGINAL_FFTW;
	else
		compres_flag = UPDATED_FFTW;
	img1->learnEpitome(compres_flag);

	////////////////////////////////////
	// Step - Reconstruction
	////////////////////////////////////
	// patchSideLengh : non overlapping extraction of xPatches
	// learn once, and reconstruct many times;
	// int newPatchSpacingNum = atoi(argv[6]); // e.g., recon-4, recon-8, different values of newPatchSpacing;
	int * i_output = new int[newPatchSpacingNum];
	img1->parseString(s_newPatchSpacing, '-', i_output);
	for (int i = 0; i != newPatchSpacingNum; ++i){
		newPatchSpacing = i_output[i];
		cout << " newPatchSpacing = " << newPatchSpacing << endl;
		for (unsigned int Read_eMean_Flag = 4; Read_eMean_Flag != 5; ++Read_eMean_Flag)
			img1->reconImgsAfterLearning(newPatchSpacing, Read_eMean_Flag);
	}


	t = clock() - t;
	printf("It took me %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);

	// ***************************************************
	// *********** release memory ************************
	// ***************************************************
	img1->clearImgEpitome();
	delete img1;
	delete[] i_output;

#endif


#define _GET_JPEGCompress_To_Epitomr_Recons_
#ifdef  _GET_JPEGCompress_To_Epitomr_Recons_

	// errorImgBit = 8, i.e., 0 - 255; errorImgBit = 4, i.e., 0 - 15;
	int errorImgBit = 8; // +9
	int jpeg_quality = 5; // +10
	string s_standard_encodeType = "JPEG2K"; // +11, e.g. = "JPEG" , or "JPEG2K"
	EncodeTypeName  standard_encodeType;
	if (s_standard_encodeType == "BMP")
		standard_encodeType = BMP;
	else if (s_standard_encodeType == "JPEG2K")
		standard_encodeType = JPEG2000;
	else if (s_standard_encodeType == "PNG")
		standard_encodeType = PNG;
	else
		standard_encodeType = JPEG;

	string encodeTypeforReconsImgs = ".bmp"; // can be changed if necessary for future usage
	string encodeTypeforErrorImgs = ".jpg"; // can be changed if necessary for future usage
	int errHistThres = 2; // +12

	QuantizeType quantizeType;
	string s_quantizetype = "NOTHING"; // +13 ;
	if (s_quantizetype == "SET2MEDIAN")
		quantizeType = SET2MEDIAN;
	else if (s_quantizetype == "MIN_MAX")
		quantizeType = BETWEEN_MIN_MAX;
	else if (s_quantizetype == "BPS_COMB4BITS_JPEG")
		quantizeType = BIT_PLANE_SLICING_COMBINE4BITSIMG_JPEG;
	else if (s_quantizetype == "BPS_JPEG")
		quantizeType = BIT_PLANE_SLICING_JPEG;
	else if (s_quantizetype == "UNI_Q_BIN")
		quantizeType = UNIFORM_QUANTIZATION_BINARY_FILE;
	else if (s_quantizetype == "UNI_Q_JPEG")
		quantizeType = UNIFORM_QUANTIZATION_JPEG;
	else{
		quantizeType = NOTHING;
		s_quantizetype = "NOTHING";
	}

	// * ReconsCompresDir="E:/imageDatabase/bmp-images-input/dog/compre-dog"
	// make the compression images dir
	MakeDir(ReconsCompresDir);
	ReconsCompresDir = ReconsCompresDir + "/" + whatkindImgs + "-EpiRecon-ErrBit-" + to_string(errorImgBit) +
		"-" + s_standard_encodeType;
	MakeDir(ReconsCompresDir);
	Distortion_Bpp * db = new Distortion_Bpp(DatabaseDir, EpitomeResultDir, ReconsCompresDir, s_standard_encodeType, standard_encodeType, quantizeType, s_quantizetype, errorImgBit, jpeg_quality,
		encodeTypeforReconsImgs, encodeTypeforErrorImgs, nameDifference, whatkindImgs, errHistThres);

	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////


	int compres_flag = 1; // +14

	// bool IsErrorImghistgram = true;// argv[mainFuncParamsIdx + 15] == "0" ? false : true; // +15
	
	// unsigned int compres_flag_2_epi_recon =  atoi(argv[mainFuncParamsIdx + 16]); // +16
	const string JPEG2K_EXE_Base_Path = "C:/Users/ccj/Desktop/CPlusPlusProjects_CCJ/openjpeg-2.1.0-win32-x86/bin"; // + 16, JPEG2K compressor and de-compressor base path
	GetJPEGCompressDir = GetJPEGCompressDir + "/" + whatkindImgs + "-" + s_standard_encodeType + "-Compress-Recons";

	// *************************************************
	// check the overhead, 
	// make sure only one parameter-setting epitome-learned result is kept in the original 
	// directory of EpitomeResultDir, i.e., EpitomeResultDir = argv[2] of the main-function;
	// *************************************************
	vector <string> EpitomeParams;
	GetDirList(EpitomeResultDir, &EpitomeParams);
	int size_EpiParams = EpitomeParams.size();
	switch (compres_flag)
	{
	case 1:
		// to get JPEG Compressed to epitome Reconstructions
		// db->getJPEGCompressToEpiRecons(IsErrorImghistgram);
		/*compres_flag_2_epi_recon = 0, 1, 2, 3*/

		/*static_cast conversion
		Converts between types using a combination of implicit and user-defined conversions.
		Syntax : static_cast < new_type > ( expression )
		*/
	{
		unsigned int compres_flag_2_epi_recon = 0;
	//	db->CompressReconToEpiRecons_plus(JPEG2K_EXE_Base_Path, static_cast <COMPRE_FLAG_2_EPI_RECON> (0));
	//	db->CompressReconToEpiRecons_plus(JPEG2K_EXE_Base_Path, static_cast <COMPRE_FLAG_2_EPI_RECON> (1));
	//	db->CompressReconToEpiRecons_plus(JPEG2K_EXE_Base_Path, static_cast <COMPRE_FLAG_2_EPI_RECON> (2));
		db->CompressReconToEpiRecons(JPEG2K_EXE_Base_Path, static_cast <COMPRE_FLAG_2_EPI_RECON> (3));
	}
		break;
	case 2:{
		// get Rate_Distortion from Directory, usually this function should be called after the above function
		// for example, we should first finish the work of JPEG-Compressed-to-epitome-Reconstructions with all compression level (JPEG: 0 - 100),
		// then, we run the function to get the rate-distortion information from the result produced above.

		// *************************************************
		// check the overhead, 
		// make sure only one parameter-setting epitome-learned result is kept in the original 
		// directory of EpitomeResultDir, i.e., EpitomeResultDir = argv[2] of the main-function;
		// *************************************************
		if (size_EpiParams == 1){
			db->getEpiCompressRate_DistortionfromDir(overhead);
		}
		else
			std::cout << "Due to the overhead, Please make sure only one parameter-setting epitome-learned result " <<
			"is kept in the original directory of EpitomeResultDir, i.e., EpitomeResultDir = argv[2] of the main-function.\n";
	}
		break;
	case 3:{
		// to get standard compression results of original input images, like JPEG, JPEG2000, etc, 
		// as the ground truth which can be used to compared with that of our method.
		// for each jpeg_quality
		db->imgCompressViaStandardCodec(GetJPEGCompressDir);
	}
		break;
	case 4:{
		// after all the possible compression level (JPEG: 0 - 100), jpeg_qualities have been finished.
		db->getStandardCodecRate_DistortionfromDir(GetJPEGCompressDir);
	}
		break;
	default:
		break;
	}



	delete db;

#endif

	return 0;
}

#endif