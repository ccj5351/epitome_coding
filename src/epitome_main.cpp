#define _epitome_fftw_convolve_time_test_
#ifdef  _epitome_fftw_convolve_time_test_
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
	string DatabaseDir = argv[1];
	string EpitomeResultDir = argv[2];
	string ReconsCompresDir = argv[3];
	string GetJPEGCompressDir = argv[3];
	bool display;  // for cout something
	if (atoi(argv[4]) != 0)
		display = true; // for cout something
	else display = false;
	string whatkindImgs = argv[5];
	const int newPatchSpacingNum = atoi(argv[6]); // e.g., recon-4, recon-8, different values of newPatchSpacing;
	const string s_newPatchSpacing = argv[7];
	// format = IsQuantizationFlag-IsDownSampleFlag-QuantizationLevel-(1/f_width)-(1/f_height)
	// f-width  : scale factor along the horizontal axis; it is computed as (double)dsize.width / src.cols;
	// f-height : scale factor along the vertical axis; it is computed as (double)dsize.height / src.rows;
	// e.g., = "0-0-100-1-1", 
	const string fileNameforRowCol_Idx_YML = argv[8]; // e.g., = "baseline"
	bool IsGrayImgs = ((atoi(argv[9])) != 0) ? true : false; // 0 : false; nonzero = true;


	// This number means the maximum images which will be saved as the final results,
	// The reason we set this parameter is that we might not save all the reconstructed images,
	// due to the disk space.
	/*the number of images will be saved as physical files*/
	int save_img_num = atoi(argv[10]); // e.g., save_img_num = 20, i.e., save 
	int nthreads = atoi(argv[11]);// thread numbers used for the parallel calculation;

	// due to the memory limit, we do parallel computation several times,
	// for each parallel computation, the OMP for(int i =0; i < MAX_PARALLEL_SIZE; ++ i)-loop 
	// can accept MAX_PARALLEL_SIZE as its ternimation limit;
	int max_omp_for_idx = atoi(argv[12]);

	// gain some important parameters from the input parameters of the main function, i.e., argv[i] ;
	int mainFuncParamsIdx = 13; // define some a number for i, just for convenience when some new parameters are added.
	int newPatchSpacing = 0;
	bool randomshifting; // make a random shift form the regular patches' positions
	string s_randomshifting = argv[mainFuncParamsIdx];
	if (stoi(s_randomshifting) != 0) // 0 : false; non-zero : true;
		randomshifting = true; // for cout something
	else randomshifting = false;
	const int epitomeWidth = atoi(argv[mainFuncParamsIdx + 1]);
	const int epitomeHeight = epitomeWidth;
	int patchSideLength = atoi(argv[mainFuncParamsIdx + 2]);
	int patchSpacing = atoi(argv[mainFuncParamsIdx + 3]);
	vector<string > inputImgsList;
	GetFileList(DatabaseDir, &inputImgsList);
	int numImgs = inputImgsList.size();
	int numIteration = atoi(argv[mainFuncParamsIdx + 4]);
	std::string s_numIteration = argv[mainFuncParamsIdx + 4];
	std::string nameDifference = argv[mainFuncParamsIdx + 5];// e.g., = "recon_"
	std::string imgEpi_EncodeType = argv[mainFuncParamsIdx + 6]; // e.g. = ".jpg", ".bmp", or ".png";
	int Epi_overheadKind_Num = atoi(argv[mainFuncParamsIdx + 7]); // file size units = Bytes, e.g., epitome Data FileSize, Units = Bytes
	// e.g., = 1, means only learn epitome;
	// e.g., = 2, means learn and reconstruction; 
	// e.g., = 3, means processing of row-column indices and direct images reconstruction;
	std::string s_flag_epi_learn_recon = argv[mainFuncParamsIdx + 8];

	// display some important input parameters of main function
	std::cout << " //*****************************************************" << endl;
	std::cout << " //***display input parameters of main function*********" << endl;
	std::cout << " //*****************************************************" << endl;
	std::cout << " Input Image Dir = " << argv[1] << endl; // #1
	std::cout << " Epitome Recons Dir = " << argv[2] << endl; // #2
	std::cout << " Recons Compression Dir = " << argv[3] << endl; // #3
	if (display != false)
		std::cout << " display = " << "true" << endl; // #4
	else
		std::cout << " display = " << "false" << endl; // #4
	std::cout << " whatkindImgs = " << whatkindImgs << endl // #5
		<< " newPatchSpacingNum = " << newPatchSpacingNum << endl // #6
		<< " s_newPatchSpacing = " << s_newPatchSpacing << endl // #7
		<< " fileNameforRowCol_Idx_YML = " << fileNameforRowCol_Idx_YML << endl; // #8
	if (IsGrayImgs) // #9
		std::cout << " IsGrayImgs = true\n";
	else
		std::cout << " IsGrayImgs = false\n";
	cout << " save_img_num = " << save_img_num << endl // #10
		<< " thread numbers for the parallel calculation = " << nthreads << endl // #11
		<< " MAX_PARALLEL_SIZE = " << max_omp_for_idx << endl;// #12


#ifdef _DEBUG
	std::cout << " recon PatchSpacing s_newPatchSpacing = " << s_newPatchSpacing << endl;
#endif // _DEBUG

	if (randomshifting != false)
		std::cout << " randomshifting = " << "true" << endl; // +0
	else
		std::cout << " randomshifting = " << "false" << endl;// +0

	std::cout << " epitomeWidth =  epitomeHeight = " << epitomeHeight << endl; // +1
	std::cout << " patchSideLengh = " << patchSideLength << endl; // +2
	std::cout << " patchSpacing = " << patchSpacing << endl; // +3
	std::cout << " EM Iteration  = " << numIteration << endl; // +4
	std::cout << " nameDifference  = " << nameDifference << endl; // +5
	std::cout << " Image Encoding Type  = " << imgEpi_EncodeType << endl; // +6
	std::cout << " Overhead Size = " << Epi_overheadKind_Num << endl; // +7
	std::cout << " Epi_learn_recon_flag = " << s_flag_epi_learn_recon << endl; // +8

	std::cout << " //*****************************************************" << endl;
	std::cout << " //*****************************************************" << endl;


#define step_reconstrcuction
#ifdef  step_reconstrcuction


	if (s_flag_epi_learn_recon == "1"){ // only learn epitome
		// ***************************************************
		// ******* initialize GrayImgEpitome Object **********
		// ***************************************************
		// make some necessary directories 
		// epitome reconstruction dir
		cout << " s_flag_epi_learn_recon == 1, that is only learn epitome, without image reconstruction.\n";
		MakeDir(EpitomeResultDir);
		EpitomeResultDir = EpitomeResultDir + "/" + whatkindImgs + "-Epi_" + s_randomshifting + "_" + std::to_string(static_cast<long long>(epitomeWidth)) + "_"
			+ std::to_string(static_cast<long long>(patchSideLength)) + "_" + std::to_string(static_cast<long long>(patchSpacing)) + "_"
			+ std::to_string(static_cast<long long>(numIteration)) + "-baseline";
		MakeDir(EpitomeResultDir);

		ImgEpitome * img1 = new ImgEpitome(
			epitomeWidth, epitomeHeight, numImgs, patchSideLength, patchSpacing, numIteration, IsGrayImgs, display,
			randomshifting, DatabaseDir, EpitomeResultDir, ReconsCompresDir, whatkindImgs, nameDifference,
			imgEpi_EncodeType, Epi_overheadKind_Num, nthreads, max_omp_for_idx);


		bool IsConstantValueInitial = false;
		double constantVal = 0.8;
		double weight_for_sigma = 0.02;
		img1->initialImgEpitome(inputImgsList, weight_for_sigma, IsConstantValueInitial, constantVal);

		////////////////////////////////////
		// Step - Learn Epitome
		////////////////////////////////////
		EpitomeVariableName compres_flag = UPDATED_FFTW;
		img1->learnEpitome(compres_flag);
		////////////////////////////////////

		// ***************************************************
		// *********** release memory ************************
		// ***************************************************
		img1->clearImgEpitome();
		delete img1;
	}

	else if (s_flag_epi_learn_recon == "2"){ // do direct baseline reconstruction;
		// but note that this case only can be successfully run when the case 1 has been run beforehand;
		// ***************************************************
		// ******* initialize GrayImgEpitome Object **********
		// ***************************************************
		// make some necessary directories 
		// epitome reconstruction dir
		cout << " s_flag_epi_learn_recon == 2, that is, to do direct baseline reconstruction, "
			<< "but note that this case will not be successfully run, "
			<< "if without the case of 's_flag_epi_learn_recon == 1' has been run beforehand.\n";
		MakeDir(EpitomeResultDir);
		EpitomeResultDir = EpitomeResultDir + "/" + whatkindImgs + "-Epi_" + s_randomshifting + "_"
			+ std::to_string(static_cast<long long>(epitomeWidth)) + "_"
			+ std::to_string(static_cast<long long>(patchSideLength)) + "_" + std::to_string(static_cast<long long>(patchSpacing)) +
			"_" + std::to_string(static_cast<long long>(numIteration)) + "-baseline";
		MakeDir(EpitomeResultDir);
		ImgEpitome * img1 = new ImgEpitome(
			epitomeWidth, epitomeHeight, numImgs, patchSideLength, patchSpacing, numIteration, IsGrayImgs, display,
			randomshifting, DatabaseDir, EpitomeResultDir, ReconsCompresDir, whatkindImgs, nameDifference,
			imgEpi_EncodeType, Epi_overheadKind_Num, nthreads, max_omp_for_idx);

		bool IsConstantValueInitial = false;
		double constantVal = 0.15;
		double weight_for_sigma = 0.01;
		img1->initialImgEpitome(inputImgsList, weight_for_sigma, IsConstantValueInitial, constantVal);

		// Step - Baseline Reconstruction
		////////////////////////////////////
		// patchSideLengh : non overlapping extraction of xPatches
		// learn once, and reconstruct many times;
		// int newPatchSpacingNum = atoi(argv[6]); // e.g., recon-4, recon-8, different values of newPatchSpacing;
		int * i_output = new int[newPatchSpacingNum];
		img1->parseString(s_newPatchSpacing, '-', i_output);
		for (int i = 0; i != newPatchSpacingNum; ++i){
			newPatchSpacing = i_output[i];
			cout << " newPatchSpacing = " << newPatchSpacing << endl;
			// here we only do Read_eMean_Flag = 1;
			unsigned int Read_eMean_Flag = 1;
			img1->reconImgsAfterLearning(newPatchSpacing, Read_eMean_Flag);
			//for (unsigned int Read_eMean_Flag = 1; Read_eMean_Flag < 2; ++Read_eMean_Flag)
			// read v_Mean from diffeerent files;
			// = 0, read v_eMean form ".yml" file;
			// = 1, read v_eMean form ".jpg" file;
			// = 2, read v_eMean form ".png" file;
			// = 3, read v_eMean form binary file;
			// = 4, read v_eMean form ".bmp" file;
			// img1->reconImgsAfterLearning(newPatchSpacing, Read_eMean_Flag);
		}
		// ***************************************************
		// *********** release memory ************************
		// ***************************************************
		img1->clearImgEpitome();
		delete img1;
		delete[] i_output;
	}

	else if (s_flag_epi_learn_recon == "3"){ // learn + baseline reconstruction
		// ***************************************************
		// ******* initialize GrayImgEpitome Object **********
		// ***************************************************
		// make some necessary directories 
		// epitome reconstruction dir
		cout << " s_flag_epi_learn_recon == 3, that is learn epitome + baseline reconstruction."
			<< "it can be considered as the combination of 's_flag_epi_learn_recon == 1' and 's_flag_epi_learn_recon == 2'.\n";
		MakeDir(EpitomeResultDir);
		EpitomeResultDir = EpitomeResultDir + "/" + whatkindImgs + "-Epi_" + s_randomshifting + "_"
			+ std::to_string(static_cast<long long>(epitomeWidth)) + "_"
			+ std::to_string(static_cast<long long>(patchSideLength)) + "_" + std::to_string(static_cast<long long>(patchSpacing)) +
			"_" + std::to_string(static_cast<long long>(numIteration)) + "-baseline";
		MakeDir(EpitomeResultDir);
		ImgEpitome * img1 = new ImgEpitome(
			epitomeWidth, epitomeHeight, numImgs, patchSideLength, patchSpacing, numIteration, IsGrayImgs, display,
			randomshifting, DatabaseDir, EpitomeResultDir, ReconsCompresDir, whatkindImgs, nameDifference,
			imgEpi_EncodeType, Epi_overheadKind_Num, nthreads, max_omp_for_idx);

		bool IsConstantValueInitial = false;
		double constantVal = 0.15;
		double weight_for_sigma = 0.01;
		img1->initialImgEpitome(inputImgsList, weight_for_sigma, IsConstantValueInitial, constantVal);


		////////////////////////////////////
		// Step - Learn Epitome
		////////////////////////////////////
		EpitomeVariableName compres_flag = UPDATED_FFTW;
		img1->learnEpitome(compres_flag);
		////////////////////////////////////

		// Step - Baseline Reconstruction
		////////////////////////////////////
		// patchSideLengh : non overlapping extraction of xPatches
		// learn once, and reconstruct many times;
		// int newPatchSpacingNum = atoi(argv[6]); // e.g., recon-4, recon-8, different values of newPatchSpacing;
		int * i_output = new int[newPatchSpacingNum];
		img1->parseString(s_newPatchSpacing, '-', i_output);
		for (int i = 0; i != newPatchSpacingNum; ++i){
			newPatchSpacing = i_output[i];
			cout << " newPatchSpacing = " << newPatchSpacing << endl;
			// here we only do Read_eMean_Flag = 1;
			unsigned int Read_eMean_Flag = 1;
			img1->reconImgsAfterLearning(newPatchSpacing, Read_eMean_Flag);
			//for (unsigned int Read_eMean_Flag = 1; Read_eMean_Flag < 2; ++Read_eMean_Flag)
			// read v_Mean from diffeerent files;
			// = 0, read v_eMean form ".yml" file;
			// = 1, read v_eMean form ".jpg" file;
			// = 2, read v_eMean form ".png" file;
			// = 3, read v_eMean form binary file;
			// = 4, read v_eMean form ".bmp" file;
			//img1->reconImgsAfterLearning(newPatchSpacing, Read_eMean_Flag);
		}
		// ***************************************************
		// *********** release memory ************************
		// ***************************************************
		img1->clearImgEpitome();
		delete img1;
		delete[] i_output;
	}

	// processing of row-column indices and get direct images reconstruction based on the baseline reconstruction;
	else if (s_flag_epi_learn_recon == "4"){

		// ***************************************************
		// ******* initialize GrayImgEpitome Object **********
		// ***************************************************
		// make some necessary directories;
		// epitome reconstruction directory;
		cout << " s_flag_epi_learn_recon == 3, that is, processing of row-column indices and direct images reconstruction,\n"
			<< " based on the baseline reconstruction.\n";
		MakeDir(EpitomeResultDir);
		string baselineEpitomeResultDir = EpitomeResultDir + "/" + whatkindImgs + "-Epi_" + s_randomshifting + "_" +
			std::to_string(static_cast<long long>(epitomeWidth)) + "_"
			+ std::to_string(static_cast<long long>(patchSideLength)) + "_" +
			std::to_string(static_cast<long long>(patchSpacing)) + "_" +
			std::to_string(static_cast<long long>(numIteration)) + "-baseline";
		EpitomeResultDir = EpitomeResultDir + "/" + whatkindImgs + "-Epi_" + s_randomshifting + "_"
			+ std::to_string(static_cast<long long>(epitomeWidth)) + "_"
			+ std::to_string(static_cast<long long>(patchSideLength)) + "_"
			+ std::to_string(static_cast<long long>(patchSpacing)) + "_"
			+ std::to_string(static_cast<long long>(numIteration)) + "-" + fileNameforRowCol_Idx_YML;
		MakeDir(EpitomeResultDir);

		ImgEpitome * img1 = new ImgEpitome(
			epitomeWidth, epitomeHeight, numImgs, patchSideLength, patchSpacing, numIteration, IsGrayImgs, display,
			randomshifting, DatabaseDir, EpitomeResultDir, ReconsCompresDir, whatkindImgs, nameDifference,
			imgEpi_EncodeType, Epi_overheadKind_Num, nthreads, max_omp_for_idx);

		////////////////////////////////////
		// Step - Reconstruction
		////////////////////////////////////

		// patchSideLengh : non overlapping extraction of xPatches
		// learn once, and reconstruct many times;
		// int newPatchSpacingNum = atoi(argv[6]); // e.g., recon-4, recon-8, different values of newPatchSpacing;

		int * i_output = new int[newPatchSpacingNum];
		img1->parseString(s_newPatchSpacing, '-', i_output);

		// format = IsQuantizationFlag-IsDownSampleFlag-QuantizationLevel-(1/f_width)-(1/f_height)
		// e.g., fileNameforRowCol_Idx_YML = "0-0-100-1-1", 
		int * i_params = new int[5];
		img1->parseString(fileNameforRowCol_Idx_YML, '-', i_params);
		bool  IsUniformQuantize = i_params[0] == 1 ? true : false; // uniform quantization flag;
		bool  IsDownSampled = i_params[1] == 1 ? true : false;  // down-sampling or up-sampling flag;
		int QuantizationLevel = i_params[2];// constant quantization level value, 0 <= QuantizationLevel <= 100;
		double f_width = 1.0 / i_params[3]; // scale factor along the horizontal axis; it is computed as (double)dsize.width / src.cols;
		double f_height = 1.0 / i_params[4];// scale factor along the vertical axis; it is computed as (double)dsize.height / src.rows;
		const bool  IsShowingImgFlag = false;
		for (int i = 0; i != newPatchSpacingNum; ++i){ // each of newPatchSpacing
			newPatchSpacing = i_output[i];
			cout << " newPatchSpacing = " << newPatchSpacing << endl;
			// here we only do Read_eMean_Flag = 1;
			unsigned int Read_eMean_Flag = 1;
			//for (unsigned int Read_eMean_Flag = 1; Read_eMean_Flag < 2; ++Read_eMean_Flag){/*read v_Mean from diffeerent files;*/
			// read v_Mean from diffeerent files;
			// = 0, read v_eMean form ".yml" file;
			// = 1, read v_eMean form ".jpg" file;
			// = 2, read v_eMean form ".png" file;
			// = 3, read v_eMean form binary file;
			// = 4, read v_eMean form ".bmp" file;

			// ***************************************************
			// do quantization and down- or up-sampling to the 
			// existing row-column indices in the format of YML files;
			// ***************************************************

			// the original row-column indices, i.e., the result from epitome learning and inference;
			string s_src_Idx = baselineEpitomeResultDir + "/"
				+ whatkindImgs + "_" +
				to_string(static_cast<long long>(newPatchSpacing)) + "_" +
				to_string(static_cast<long long>(Read_eMean_Flag)) + "_RowCol_Idx.yml";
			// the quantized indices of the original ones;
			string s_quantized_Idx = EpitomeResultDir + "/" + "quantized-RowCol-Idx";
			MakeDir(s_quantized_Idx);
			s_quantized_Idx += "/" + whatkindImgs + "_" +
				to_string(static_cast<long long>(newPatchSpacing)) + "_"
				+ to_string(static_cast<long long>(Read_eMean_Flag)) + "-quanti_RowCol_Idx.yml";
			// the reconstruction of the quantized indices via de-quantization and/or down- or up-sampling.
			string s_dst_Idx = EpitomeResultDir + "/"
				+ whatkindImgs + "_" +
				to_string(static_cast<long long>(newPatchSpacing)) + "_"
				+ to_string(static_cast<long long>(Read_eMean_Flag)) + "_RowCol_Idx.yml";


			std::cout << "do quantization and down- or up-sampling to the existing row-column indices in the format of YML files\n"
				<< " QuantizationLevel = " << QuantizationLevel << endl
				<< " Calculating...\n";
			UniformQuant_row_col_idx_saved_as_YML(whatkindImgs, s_src_Idx, s_dst_Idx, s_quantized_Idx, patchSideLength, newPatchSpacing, epitomeWidth,
				IsShowingImgFlag, IsUniformQuantize, (double)QuantizationLevel, IsDownSampled, f_width, f_height);
			std::cout << "Finish it.\n";
			std::cout << "Start doing reconImgsDirect(...).\n";
			img1->reconImgsDirect(newPatchSpacing, baselineEpitomeResultDir, Read_eMean_Flag);
			std::cout << "Finish doing reconImgsDirect(...).\n";
			//}/*end of reading v_Mean from diffeerent files;*/
		} /*end of each of newPatchSpacing*/

		// ***************************************************
		// *********** release memory ************************
		// ***************************************************
		img1->clearImgEpitome();
		delete img1;
		delete[] i_output;
		delete[] i_params;

	} /*end of s_flag_epi_learn_recon == "3" */

	else {
		std::cout << "Do Nothing! s_flag_epi_learn_recon should be 1, 2 or 3.\n";
	}

	t = clock() - t;
	printf("It took me %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);



#endif


#define _GET_Compress_To_Epitomr_Recons_
#ifdef  _GET_Compress_To_Epitomr_Recons_

	// errorImgBit = 8, i.e., 0 - 255; errorImgBit = 4, i.e., 0 - 15;
	int errorImgBit = atoi(argv[mainFuncParamsIdx + 9]); // +9
	int jpeg_quality = atoi(argv[mainFuncParamsIdx + 10]); // +10
	string s_standard_encodeType = argv[mainFuncParamsIdx + 11]; // +11, e.g. = "JPEG" , or "JPEG2K"

	string encodeTypeforReconsImgs = ".bmp"; // can be changed if necessary for future usage

	string encodeTypeforErrorImg; // can be changed if necessary for future usage

	EncodeTypeName  standard_encodeType;
	if (s_standard_encodeType == "BMP"){
		standard_encodeType = BMP;
		encodeTypeforErrorImg = ".jpg";
	}
	else if (s_standard_encodeType == "JPEG2K"){
		standard_encodeType = JPEG2000;
		encodeTypeforErrorImg = ".j2k";
	}
	else if (s_standard_encodeType == "PNG"){
		standard_encodeType = PNG;
		encodeTypeforErrorImg = ".png";
	}

	else{
		standard_encodeType = JPEG;
		encodeTypeforErrorImg = ".jpg";
	}


	int errHistThres = atoi(argv[mainFuncParamsIdx + 12]); // +12

	QuantizeType quantizeType;
	string s_quantizetype = argv[mainFuncParamsIdx + 13]; // +13 ;
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
	ReconsCompresDir = ReconsCompresDir + "/" + whatkindImgs + "-EpiRecon-ErrBit-"
		+ to_string(static_cast<long long>(errorImgBit)) +
		"-" + s_standard_encodeType + "-" + fileNameforRowCol_Idx_YML;
	MakeDir(ReconsCompresDir);
	Distortion_Bpp * db = new Distortion_Bpp(DatabaseDir, EpitomeResultDir, ReconsCompresDir, s_standard_encodeType, standard_encodeType, quantizeType, s_quantizetype, errorImgBit, jpeg_quality,
		encodeTypeforReconsImgs, encodeTypeforErrorImg, nameDifference, whatkindImgs, errHistThres, save_img_num, IsGrayImgs);
	db->setImgReadFlag();
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////


	int compres_flag = atoi(argv[mainFuncParamsIdx + 14]); // +14
	// bool IsErrorImghistgram = argv[mainFuncParamsIdx + 15] == "0" ? false : true; // +15, currently this parameter is not used.
	const string JPEG2K_EXE_Base_Path = argv[mainFuncParamsIdx + 15]; // + 15, JPEG2K compressor and de-compressor base path
	GetJPEGCompressDir = GetJPEGCompressDir + "/" + whatkindImgs + "-" + s_standard_encodeType + "-Direct-Compress-Recon";

	// *************************************************
	// check the overhead, 
	// make sure only one parameter-setting epitome-learned result is kept in the original 
	// directory of EpitomeResultDir, i.e., EpitomeResultDir = argv[2] of the main-function;
	// *************************************************

	switch (compres_flag){

	case 1: // To calculate the overhead in order to calculate the rate-distortion result,
		// specifically, to change row_column indices into j2k image and/or ".7z" files for a less compressed file size,
		// and the compressed indices is considered as the overhead;
		// the overhead is a necessary parameter for the following code "case 3" to get Rate_Distortion.
	{
		// e.g., C:/Users/ccj/imageDatabase/bmp-images-epi-result/epi-6_Try_epi-sub-Edwin-less-full/1-0-030-2-1;
		const string basePath_Idx = argv[mainFuncParamsIdx + 17];
		const string J2K_Result_BasePath = basePath_Idx;
		vector<string> RowCol_IdxFilelist;
		GetFileList(basePath_Idx, &RowCol_IdxFilelist);

		string fileStreamForIdxSize = J2K_Result_BasePath + "/" + "Idx-J2K-Size";
		MakeDir(fileStreamForIdxSize);
		fileStreamForIdxSize += "/" + whatkindImgs + "-Idx-J2K-Size.txt";
		ofstream fout_IdxImgSize(fileStreamForIdxSize, std::ofstream::app);
		if (!fout_IdxImgSize){
			std::cout << "File Not Opened" << endl;
		}
		else{
			auto SizeIdxFilelist = RowCol_IdxFilelist.size();
			vector <uintmax_t> v_size_7z_file(SizeIdxFilelist),
				v_compreRowColIdx_fileSize(SizeIdxFilelist);

			for (int i = 0; i != SizeIdxFilelist; ++i){
				std::size_t pos = RowCol_IdxFilelist[i].find(".");
				string category_RowColIdx = RowCol_IdxFilelist[i].substr(0, pos); // e.g., "Edwin-less_4_0_RowCol_Idx";
				string path_Idx = basePath_Idx + "/" + RowCol_IdxFilelist[i];
				const int patchSpace = atoi(category_RowColIdx.substr(whatkindImgs.length() + 1, 1).c_str());
				v_compreRowColIdx_fileSize[i] =
					saveRowColIdx2_J2KImg(patchSideLength, patchSpace, epitomeWidth, whatkindImgs, JPEG2K_EXE_Base_Path, category_RowColIdx, path_Idx,
					J2K_Result_BasePath, true, v_size_7z_file[i]);
				fout_IdxImgSize << category_RowColIdx << " J2K images size = " << v_compreRowColIdx_fileSize[i]
					<< " Bytes, 7z file size = " << v_size_7z_file[i] << " Bytes;\n";
			}
			// for using the data conveniently,
			// at the end, print them together in one line,
			// separated by comma;
			fout_IdxImgSize << "J2K row-column indices file sizes (Bytes) are \n";
			for (int i = 0; i != SizeIdxFilelist; ++i){
				fout_IdxImgSize << v_compreRowColIdx_fileSize[i] << ",";
			}
			fout_IdxImgSize << "\n7z row-column indices file size (Bytes) are \n";
			for (int i = 0; i != SizeIdxFilelist; ++i){
				fout_IdxImgSize << v_size_7z_file[i] << ",";
			}
			fout_IdxImgSize << endl;
		}
		fout_IdxImgSize.close();
		vector<string>().swap(RowCol_IdxFilelist);
	}
	break;

	case 2:
		// Get the final reconstruction result of our method,
		// i.e., to do compression to the residual between the epitome reconstructions and the original input images;
		/*compres_flag_2_epi_recon = 0, 1, 2, 3*/
	{
		// Previous function, currently they are not used. 
		// bool IsErrorImghistgram = argv[mainFuncParamsIdx + 15] == "0" ? false : true; // +15
		// db->getJPEGCompressToEpiRecons(IsErrorImghistgram);
		unsigned int compres_flag_2_epi_recon = atoi(argv[mainFuncParamsIdx + 16]); // should be 0, 1, 2, 3;
		// compres_flag_2_epi_recon  = 
		// 0, means UNIFORM_QUANTIZATION;
		// 1, means UNIFORM_QUANTIZATION_LOSSLESS_J2K;
		// 2, means UNIFORM_QUANTIZATION_LOSSY_JPEG;
		// 3, means UNIFORM_QUANTIZATION_LOSSY_J2K;

		/*static_cast conversion
		Converts between types using a combination of implicit and user-defined conversions.
		Syntax : static_cast < new_type > ( expression )
		*/
		db->CompressReconToEpiRecons(JPEG2K_EXE_Base_Path, static_cast <COMPRE_FLAG_2_EPI_RECON> (compres_flag_2_epi_recon));
	}
	break;
	case 3:{
		// get Rate_Distortion from Directory, usually this function should be called after the above function
		// for example, we should first finish the work of Compress-to-epitome-Reconstructions with all compression level (e.g., JPEG: 0 - 100),
		// then, we run this function to get the rate-distortion curves from the result produced above.

		// *************************************************
		// check the overhead, 
		// make sure only one parameter-setting epitome-learned result is kept in the original 
		// directory of EpitomeResultDir, i.e., EpitomeResultDir = argv[2] of the main-function;
		// *************************************************

		const uintmax_t sum_input_size = db->getFileSizeFromDir(DatabaseDir);
		std::map <string, uintmax_t> map_Epi_overhead;
		vector<string> v_key(Epi_overheadKind_Num);

		// i-j: patchSpacing & read_eMean_flag;
		// patchSpacing = 4, 5, 6, 7, 8;
		// read_eMean_flag = 0, 1, 2, 3, 4(maybe not included sometimes);
		// e.g., 4-0, 4-1, 4-2, 4-3;
		// e.g., 8-0, 8-1, 8-2, 8-3;
		for (int i = 4, idx = 0; i != 9; ++i)
			for (int j = 0; j != 4; ++j, ++idx)
				v_key[idx] =
				to_string(static_cast<long long>(i)) + "_" +
				to_string(static_cast<long long>(j));
		// *************************************************
		// epitome indices file size for 20 categories
		// they should be changed in different experiment cases;
		// unit is Bytes;
		uintmax_t eMeanFileSize_0 = 79289, // 0 = compressed ".yml" files, e.g., ".yml" --> ".7z"
			eMeanFileSize_1 = 5987, // 1 = ".jpg"
			eMeanFileSize_2 = 11636, // 2 = ".png"
			eMeanFileSize_3 = 79289; // 3 = compressed ".bin" files, e.g., ".bin" --> ".7z"

		int * v_Epi_Idx = new int[Epi_overheadKind_Num];
		string s_Epi_Idx_FileSize = argv[mainFuncParamsIdx + 18];
		db->parseString(s_Epi_Idx_FileSize, ',', v_Epi_Idx);

		uintmax_t Epi_overheads[] = {
			v_Epi_Idx[0] + eMeanFileSize_0, v_Epi_Idx[1] + eMeanFileSize_1, v_Epi_Idx[2] + eMeanFileSize_2, v_Epi_Idx[3] + eMeanFileSize_3,
			v_Epi_Idx[4] + eMeanFileSize_0, v_Epi_Idx[5] + eMeanFileSize_1, v_Epi_Idx[6] + eMeanFileSize_2, v_Epi_Idx[7] + eMeanFileSize_3,
			v_Epi_Idx[8] + eMeanFileSize_0, v_Epi_Idx[9] + eMeanFileSize_1, v_Epi_Idx[10] + eMeanFileSize_2, v_Epi_Idx[11] + eMeanFileSize_3,
			v_Epi_Idx[12] + eMeanFileSize_0, v_Epi_Idx[13] + eMeanFileSize_1, v_Epi_Idx[14] + eMeanFileSize_2, v_Epi_Idx[15] + eMeanFileSize_3,
			v_Epi_Idx[16] + eMeanFileSize_0, v_Epi_Idx[17] + eMeanFileSize_1, v_Epi_Idx[18] + eMeanFileSize_2, v_Epi_Idx[19] + eMeanFileSize_3
		};

		// C++ 98
		vector <uintmax_t> v_Epi_overhead(Epi_overheads, Epi_overheads + sizeof(Epi_overheads) / sizeof(uintmax_t));
		// *************************************************
		for (int k = 0; k != Epi_overheadKind_Num; ++k)
			map_Epi_overhead[v_key[k]] = v_Epi_overhead[k];
		db->getEpiCompressRate_DistortionfromDir(map_Epi_overhead, sum_input_size);
		// std::cout << "Due to the overhead, Please make sure only one parameter-setting epitome-learned result " <<
		// "is kept in the original directory of EpitomeResultDir, i.e., EpitomeResultDir = argv[2] of the main-function.\n";
	}
		   break;
	case 4:{
		// to get standard compression results of original input images, like JPEG, JPEG2000, etc, 
		// as the ground truth which can be used to compared with that of our method.
		// for each jpeg_quality
		db->imgCompressViaStandardCodec(GetJPEGCompressDir, JPEG2K_EXE_Base_Path);
	}
		   break;
	case 5:{
		// after all the possible compression level (JPEG: 0 - 100), jpeg_qualities have been finished.
		db->getStandardCodecRate_DistortionfromDir(GetJPEGCompressDir);
	}
		   break;
	default:
		std::cout << "Do nothing for compres_flag!\n";
		break;
	}

	delete db;

#endif

	return 0;
}

#endif