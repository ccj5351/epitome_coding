#include "ImageEpitome.hpp"

// constructor
ImgEpitome::ImgEpitome(int width, int height, int numImgs, int patchLengh, int patchSpac, int numItera,
	bool isGray, bool display, bool random,
	std::string DatabaseDir,
	std::string EpitomeResultDir,
	std::string ReconsCompresDir,
	std::string whatkindImgs,
	std::string nameDifference,// e.g., = = "recon_"
	std::string imgEncodeType, // e.g. = ".jpg", ".bmp", or ".png";
	int overhead,// file size of epitome overhead, Units = Bytes
	int nthreads, // thread numbers used for the parallel calculation;
	int max_omp_for_idx
	)  
	: eWidth(width),
	eHeight(height),
	patchSideLengh(patchLengh),
	patchSpacing(std::max<int>(patchSpac, 2)), // the distance between sampled patches from the input image.  it should be >= 2
	numIteration(numItera),
	isgrayScale(isGray),
	verbose(display),
	RandomShifting(random),
	DatabaseDir(DatabaseDir),
	EpitomeResultDir(EpitomeResultDir),
	ReconsCompresDir(ReconsCompresDir),
	whatkindImgs(whatkindImgs),
	nameDifference(nameDifference),
	imgEncodeType(imgEncodeType),
	overhead(overhead)
{
	//this->v_eMean = new double[eWidth*eHeight]; //the Epitome, eMean, all zeros as initialized values
	//this->v_eVar = new double[eWidth*eHeight]; //the Epitome, eVar , all ones as initialized values
	this->v_maxPost_row_column_Idx = vector<vector<int>>(numImgs, vector<int>()); // indices learned from EM algorithm and will be saved into files (e.g.".yml")
	//	this->v_row_Idx = vector<vector<int>>(numImgs, vector<int>()); // row indices of xPatches
	//	this->v_col_Idx = vector<vector<int>>(numImgs, vector<int>()); // column indices of xPatches
	// initialize the variance to unity;
	int eNum = width * height;
	if (isgrayScale){
		this -> Epi_Image_Read_Flag = 0;
		// initialize the variance to ones;
		this->v_eVar_Gray = vector<Epitome_Doulbe>(eNum, 1);
		// initialize the mean to zeros;
		this->v_eMean_Gray = vector<Epitome_Doulbe>(eNum, 0);
	}
	else{
		this -> Epi_Image_Read_Flag = 1;
		this->v_eVar_Red = vector<Epitome_Doulbe>(eNum, 1);
		this->v_eVar_Gre = vector<Epitome_Doulbe>(eNum, 1);
		this->v_eVar_Blu = vector<Epitome_Doulbe>(eNum, 1);
		// initialize the mean to zeros;
		this->v_eMean_Red = vector<Epitome_Doulbe>(eNum, 0);
		this->v_eMean_Gre = vector<Epitome_Doulbe>(eNum, 0);
		this->v_eMean_Blu = vector<Epitome_Doulbe>(eNum, 0);
	}

	this->NUM_THREADS = nthreads;
	this->MAX_PARALLEL_SIZE = max_omp_for_idx;
	v_double_sumQ = vector<Epitome_Doulbe>(eNum, 0.);
	if (isgrayScale){
		v_double_sumQX_Gray = vector<Epitome_Doulbe>(eNum, 0.);
		v_double_sumQXX_Gray = vector<Epitome_Doulbe>(eNum, 0.);
	}
	else{
		v_double_sumQX_R = vector<Epitome_Doulbe>(eNum, 0.);
		v_double_sumQXX_R = vector<Epitome_Doulbe>(eNum, 0.);
		v_double_sumQX_G = vector<Epitome_Doulbe>(eNum, 0.);
		v_double_sumQXX_G = vector<Epitome_Doulbe>(eNum, 0.);
		v_double_sumQX_B = vector<Epitome_Doulbe>(eNum, 0.);
		v_double_sumQXX_B = vector<Epitome_Doulbe>(eNum, 0.);
	}
}



void ImgEpitome::displayEpitome()
{
	std::cout << "Information of Epitome::\n";
	std::cout << "the size of Epitome is : " << eHeight << " by " << eWidth << endl
		<< "patchSideLengh = " << patchSideLengh << " ,patchSpacing = " << patchSpacing << endl;

	if (verbose != false)
		std::cout << "displayProcess = " << "true";
	else
		std::cout << "displayProcess = " << "false";

	if (RandomShifting != false)
		std::cout << ", RandomShifting = " << "true" << endl;
	else
		std::cout << " , RandomShifting = " << "false" << endl;

	std::cout << "numIteration = " << numIteration << endl
		<< "eMean has been randomized in the interval [0, 1], and eVar has been set to ones.\n";
}


void ImgEpitome::updateEpitome(const char * value, EpitomeVariableName type){
	switch (type)
	{
	case IMG_ENCODE_TYPE:
		imgEncodeType = value;
		break;
	case NAME_DIFFERENCE:
		nameDifference = value;
		break;
	case OVERHEAD_SIZE:
		overhead = atoi(value);
		break;
	default:
		cout << "Parameter type does not exist, please check it!\n";
		break;
	}
}



// normalize array, from [-a, a] to [2*newMedian-newMaxVal, newMaxVal]
void ImgEpitome::normalize(Epitome_Doulbe * input, int num, Epitome_Doulbe minVal, Epitome_Doulbe maxVal, Epitome_Doulbe newMedian, Epitome_Doulbe newMaxVal){
	Epitome_Doulbe absMax = max<Epitome_Doulbe>(abs(minVal), abs(maxVal));

	for (int i = 0; i != num; ++i){
		input[i] = (input[i] / absMax) * (newMaxVal - newMedian) + newMedian;
	}

	// pay attention the initialization value here
	double minCoeff = DOUBLE_MAX;
	double maxCoeff = DOUBLE_MIN;

	for (int i = 0; i != num; ++i){
		minCoeff = minCoeff < input[i] ? minCoeff : input[i];
		maxCoeff = maxCoeff > input[i] ? maxCoeff : input[i];
	}
	if (ImgEpitome::verbose){
		if (maxCoeff == newMaxVal){
			std::cout << "\nMinCoeff = " << minCoeff << "  , MaxCoeff = " << maxCoeff << endl;
			std::cout << "Normalization has been done, which means that all the elements are in the closed interval of [" << 2 * newMedian - newMaxVal << " , "
				<< newMaxVal << " ].\n";
		}
		else
			std::cout << "Normalization has NOT been done!!\n" << endl;
	}
}


// initialized Epitome Mean
// the means  for each color channel in the epitome, were initialized to the mean of all values in the same channel
// in the training set, plus Gaussian noise with 1/100th of the standard deviation (i.e., sigma) of the training data.
// The variances were initialized to the variance of the training data.
void ImgEpitome::randomEpitomeMean(vector<string> & v_ImgNames, // just image names (including file extension), without full directory in this parameter
	Epitome_Doulbe weight_for_sigma // e.g., = 1/100 
	// (see "plus Gaussian noise with 1/100th of the standard deviation (i.e., sigma) of the training data")
	){
	// in OpenCV, the 0-channel is Blue, 1-channel is Green, and the last one is the Red.
	// the means  for each color channel in the epitome, were initialized to the mean of all values in the same channel
	// in the training set, plus Gaussian noise with 1/100th of the standard deviation (i.e., sigma) of the training data.
	// firstly, we want to calculate a mean and variance based on the current image, 

	// to calculate pixelMean and pixelStd
	Epitome_Doulbe sumX_R = 0.,  // sum of all the pixels
		sumXX_R = 0., // sum of the squares of pixels
		sumX_G = 0.,  // sum of all the pixels
		sumXX_G = 0., // sum of the squares of pixels
		sumX_B = 0.,  // sum of all the pixels
		sumXX_B = 0., // sum of the squares of pixels
		sumX_Gray = 0.,  // sum of all the pixels
		sumXX_Gray = 0.; // sum of the squares of pixels
	long int sumPixel = 0;

	for (int imgIdx = 0, idxSize = v_ImgNames.size(); imgIdx != idxSize; ++imgIdx){
		//C++: Mat::Mat(const Mat& m),opencv Mat constructors;
		// parameter: m ¨C Array that (as a whole or partly) is assigned to the constructed matrix. 
		// No data is copied by these constructors. Instead, the header pointing to m data or its sub-array is constructed and associated with it.
		// The reference counter, if any, is incremented. 
		// So, when you modify the matrix formed using such a constructor, you also modify the corresponding elements of m . 
		// If you want to have an independent copy of the sub-array, use Mat::clone().
		// Copies the matrix to another one.
		// C++: void Mat::copyTo(OutputArray m) const
		cv::Mat_<double> Img = (cv::imread(DatabaseDir + "/" + v_ImgNames[imgIdx], Epi_Image_Read_Flag)) / 255.0;
		// size of input image
		sumPixel += (Img.cols) * (Img.rows);
		cv::Mat_<double> Img2 = Img.mul(Img); // I^2
		if (Img.empty()){
			cerr << " > No input image is being read!" << endl;
		}
		else{
			// Template class for a 4-element vector derived from Vec.
			// typedef Scalar_<double> Scalar;
			// for gray images: Scalar::channel[0] is valid; channel[1] ~[3] = 0;
			// for color images: Scalar::channel[0]~[2] is valid; channel[3] = 0;
			cv::Scalar s = cv::sum(Img);         // sum elements per channel
			cv::Scalar s2 = cv::sum(Img2);         // sum elements per channel
			if (isgrayScale){ /*gray images*/
				sumX_Gray += s.val[0]; //
				sumXX_Gray += s2.val[0]; //
			}
			else{ /*color images*/
				sumX_R += s.val[2]; // red channel
				sumX_G += s.val[1]; // green channel
				sumX_B += s.val[0]; // blue channel
				sumXX_R += s2.val[2]; // red channel
				sumXX_G += s2.val[1]; // green channel
				sumXX_B += s2.val[0]; // blue channel
			}

		} /*if image is successfully read.*/
	} /*end of each image*/
	if (!isgrayScale){ /*color images*/
		std::cout << " sumX_R = " << sumX_R << " , sumX_G = " << sumX_G << " , sumX_B = " << sumX_B << endl;
		std::cout << " sumXX_R = " << sumXX_R << " , sumXX_G = " << sumXX_G << " , sumXX_B = " << sumXX_B << endl;
	}
	else{/*gray images*/
		std::cout << " sumX_Gray = " << sumX_Gray << " , and sumXX_Gray = " << sumXX_Gray << endl;
	}

	// mean = (sum) divided by (# of the elements)
	// 0: Red; 1: Green; 2: Blue; 3: Gray-scale;
	Epitome_Doulbe pixelMean[4] = { sumX_R / sumPixel, sumX_G / sumPixel, sumX_B / sumPixel, sumX_Gray / sumPixel };
	// variance = 1/n*sum((xi - mean)^2) = 1/n*sum(xi^2) - mean^2
	Epitome_Doulbe pixelStd[4] = {
		sqrt(sumXX_R / sumPixel - pixelMean[0] * pixelMean[0]),
		sqrt(sumXX_G / sumPixel - pixelMean[1] * pixelMean[1]),
		sqrt(sumXX_B / sumPixel - pixelMean[2] * pixelMean[2]),
		sqrt(sumXX_Gray / sumPixel - pixelMean[3] * pixelMean[3]), };
	std::cout << " pixelMean(R, G, B, Gray) = " << pixelMean[0] << ", " << pixelMean[1] << ", " << pixelMean[2] << ", " << pixelMean[3]
		<< " pixelStd(R, G, B, Gray) = " << pixelStd[0] << ", " << pixelStd[1] << ", " << pixelStd[2] << ", " << pixelStd[3] << endl;

	// pay attention the initialization value here
	Epitome_Doulbe minCoeff[4] = { DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX, DOUBLE_MAX };
	Epitome_Doulbe maxCoeff[4] = { DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN, DOUBLE_MIN };

	// make the initialized mean normally distributed using the sample mean,
	// and standard deviation;
	// but making sure to keep the pixels between 0 and 1
	if (isgrayScale){ /*gray image*/
		// get Normal_distribution to v_eMean_Gray;
		getNormal_distribution(v_eMean_Gray, pixelMean[3], pixelStd[3], minCoeff[3], maxCoeff[3], weight_for_sigma);
		// if eMean exceeds [0, 1], then do normalization of Mean of the Epitome
		if ((minCoeff[3] < 0) | (maxCoeff[3] > 1))
			normalize(&v_eMean_Gray[0], eHeight*eWidth, minCoeff[3], maxCoeff[3], 0.5, 1.0);
	}
	else {/*color image*/
		// get Normal_distribution to v_eMean_Red,v_eMean_Gre, and v_eMean_Blu ;
		getNormal_distribution(v_eMean_Red, pixelMean[0], pixelStd[0], minCoeff[0], maxCoeff[0], weight_for_sigma);
		getNormal_distribution(v_eMean_Gre, pixelMean[1], pixelStd[1], minCoeff[1], maxCoeff[1], weight_for_sigma);
		getNormal_distribution(v_eMean_Blu, pixelMean[2], pixelStd[2], minCoeff[2], maxCoeff[2], weight_for_sigma);
		// if eMean exceeds [0, 1], then do normalization of Mean of the Epitome
		if ((minCoeff[2] < 0) | (maxCoeff[2] > 1))
			normalize(&v_eMean_Blu[0], eHeight*eWidth, minCoeff[2], maxCoeff[2], 0.5, 1.0);
		if ((minCoeff[1] < 0) | (maxCoeff[1] > 1))
			normalize(&v_eMean_Gre[0], eHeight*eWidth, minCoeff[1], maxCoeff[1], 0.5, 1.0);
		if ((minCoeff[0] < 0) | (maxCoeff[0] > 1))
			normalize(&v_eMean_Red[0], eHeight*eWidth, minCoeff[0], maxCoeff[0], 0.5, 1.0);
	}

	cout << "Normal_distribution of all the eMeans (i.e., R, G, B and Gray) finishes.\n";
}

void ImgEpitome::initialImgEpitome(vector<string> & v_ImgPath, // image names, including full directory in this parameter
	double weight_for_sigma, // e.g., = 1/100 
	// (see "plus Gaussian noise with 1/100th of the standard deviation (i.e., sigma) of the training data")
	bool IsConstantValueInitial, 
	double constantVal
	){
	int i_size = eHeight * eWidth;
	// randomize the eMean;
	if (!IsConstantValueInitial){
		randomEpitomeMean(v_ImgPath, weight_for_sigma);
	}
	else{ /*initialization as constants */
		if (isgrayScale){
			for (int i = 0; i != i_size; ++i)
				v_eMean_Gray[i] = constantVal;
			// initial the eVar as ones
			for (int i = 0; i != i_size; ++i)
				v_eVar_Gray[i] = 1.0;
		}
		else{
			for (int i = 0; i != i_size; ++i)
				v_eMean_Red[i] = constantVal;
			for (int i = 0; i != i_size; ++i)
				v_eMean_Gre[i] = constantVal;
			for (int i = 0; i != i_size; ++i)
				v_eMean_Blu[i] = constantVal;
			// initial the eVar as ones
			for (int i = 0; i != i_size; ++i)
				v_eVar_Red[i] = 1.0;
			for (int i = 0; i != i_size; ++i)
				v_eVar_Gre[i] = 1.0;
			for (int i = 0; i != i_size; ++i)
				v_eVar_Blu[i] = 1.0;
		}
		std::cout << "eMeans has been set to all the values of " << constantVal << ";" << endl;
		std::cout << "eVars has been set to all the values of 1.0;" << endl;
	}
	ImgEpitome::displayEpitome();
}


void ImgEpitome::clearImgEpitome(){
	
	// release the eMean and eVar
	if (isgrayScale){/*gray images*/
		// initialize the variance to ones;
		vector<Epitome_Doulbe>().swap(v_eVar_Gray);
		vector<Epitome_Doulbe>().swap(v_eMean_Gray);
	}
	else{/*color images*/
		vector<Epitome_Doulbe>().swap(v_eMean_Red);
		vector<Epitome_Doulbe>().swap(v_eVar_Red);
		vector<Epitome_Doulbe>().swap(v_eMean_Gre);
		vector<Epitome_Doulbe>().swap(v_eVar_Gre);
		vector<Epitome_Doulbe>().swap(v_eMean_Blu);
		vector<Epitome_Doulbe>().swap(v_eVar_Blu);
	}

	// release 
	vector<Epitome_Doulbe>().swap(v_double_sumQ);
	vector<Epitome_Doulbe>().swap(v_double_sumQX_R);
	vector<Epitome_Doulbe>().swap(v_double_sumQXX_R);
	vector<Epitome_Doulbe>().swap(v_double_sumQX_G);
	vector<Epitome_Doulbe>().swap(v_double_sumQXX_G);
	vector<Epitome_Doulbe>().swap(v_double_sumQX_B);
	vector<Epitome_Doulbe>().swap(v_double_sumQXX_B);
	vector<Epitome_Doulbe>().swap(v_double_sumQX_Gray);
	vector<Epitome_Doulbe>().swap(v_double_sumQXX_Gray);

	vector<Epitome_Doulbe>().swap(v_patchSize_Ones);
	/*R channel*/
	vector<Epitome_Doulbe>().swap(v_eMeanOverVar_R);
	vector<Epitome_Doulbe>().swap(v_InvVar_R);
	vector<Epitome_Doulbe>().swap(v_LogVar_R); // not xPatch
	vector<Epitome_Doulbe>().swap(v_eMean2OverVar_R); // not xPatch
	vector<Epitome_Doulbe>().swap(v_eLogVarSum_R);
	vector<Epitome_Doulbe>().swap(v_eeSum_R);
	
	/*G channel*/
	vector<Epitome_Doulbe>().swap(v_eMeanOverVar_G);
	vector<Epitome_Doulbe>().swap(v_InvVar_G);
	vector<Epitome_Doulbe>().swap(v_LogVar_G);// not xPatch
	vector<Epitome_Doulbe>().swap(v_eMean2OverVar_G); // not xPatch
	vector<Epitome_Doulbe>().swap(v_eLogVarSum_G);
	vector<Epitome_Doulbe>().swap(v_eeSum_G);
	/*B channel*/
	vector<Epitome_Doulbe>().swap(v_eMeanOverVar_B);
	vector<Epitome_Doulbe>().swap(v_InvVar_B);
	vector<Epitome_Doulbe>().swap(v_LogVar_B); // not xPatch
	vector<Epitome_Doulbe>().swap(v_eMean2OverVar_B); // not xPatch
	vector<Epitome_Doulbe>().swap(v_eLogVarSum_B);
	vector<Epitome_Doulbe>().swap(v_eeSum_B);
	/*Gray channel*/
	vector<Epitome_Doulbe>().swap(v_eMeanOverVar_Gray);
	vector<Epitome_Doulbe>().swap(v_InvVar_Gray);
	vector<Epitome_Doulbe>().swap(v_LogVar_Gray); // not xPatch
	vector<Epitome_Doulbe>().swap(v_eMean2OverVar_Gray); // not xPatch
	vector<Epitome_Doulbe>().swap(v_eLogVarSum_Gray);
	vector<Epitome_Doulbe>().swap(v_eeSum_Gray);

	// release the indices
	vector<vector<int>>().swap(v_maxPost_row_column_Idx);

}


void ImgEpitome::parseString(const string & s_input, // e.g., == "4-8-12"
	const char & flag, // e.g., = '-'
	int * i_output // output which have been parsed based on the input string
	){

	// function -- std::string::length, returns the length of the string, in terms of bytes.
	string s_temp = s_input;
	const int length = s_input.length();
	int i = 0;
	std::size_t pos = 0;
	while (s_temp.find(flag) != std::string::npos){ // while there exists "flag"

		// do the loop, when we can find the 'flag'
		// function -- std::basic_string::find;
		// Returns position of the first character of the found substring or "npos" if no such substring is found.
		// std::basic_string::npos
		cout << " preprocessed s_temp = " << s_temp << endl;
		pos = s_temp.find(flag);      // e.g., position of ".jpg" in nameTemp
		*(i_output + i) = stoi(s_temp.substr(0, pos));     // get beginning to ".jpg"
		s_temp = s_temp.substr(pos + 1, length + 1);
		cout << " processed s_temp = " << s_temp << endl;
		++i;
	}
	// for the last charts after the lase flag, like, flag = '-'
	*(i_output + i) = stoi(s_temp);
#ifdef _DEBUG
	std::cout << " i_output = ";
	for (int j = 0; j <= i; ++j){
		std::cout <<  *(i_output + j) << " ";
	}
	std::cout << endl;
#endif // _DEBUG 


}

vector<vector<int>>  ImgEpitome::read_row_col_idx_form_YML(
	const string & s_Idx //  e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
	){

	FileStorage fs_Idx_out(s_Idx, FileStorage::READ);
	vector<vector<int>> v_maxPost_row_column_Idx;
	if (!fs_Idx_out.isOpened()){
		cerr << "failed to open " << s_Idx << endl;
	}
	else{
		FileNode fn_Idx = fs_Idx_out[whatkindImgs];
		FileNodeIterator it = fn_Idx.begin(), it_end = fn_Idx.end();
		int numInputImage = fn_Idx.size(); // the number of images
		int temp_idx = 0;
		v_maxPost_row_column_Idx = vector<vector<int>>(numInputImage, vector<int>());
		for (; it != it_end; ++it, temp_idx++){
			// v_maxPost_row_column_Idx[temp_idx].clear(); // max-posterior indices read from the beforehand saved files (e.g.".yml")
			(*it)["maxPost_row_col_Idx"] >> v_maxPost_row_column_Idx[temp_idx];
		}
	}
	fs_Idx_out.release();
	return v_maxPost_row_column_Idx;

}

void  ImgEpitome::write_row_col_idx_2_YML(
	const string & s_dst_Idx, //  e.g., =  "E:/imageDatabase/bmp-images-input/dog/row-col-idx/Edwin-less_4_RowCol_Idx.yml";
	int & numInputImage // the number of images
	){
	//******************************************************************
	// third, save the uniform quantized and/or down-sampled/up-sampled 
	// indices to the output YML file;
	//******************************************************************
	FileStorage fs_Idx_write(s_dst_Idx, FileStorage::WRITE);
	if (!fs_Idx_write.isOpened()){
		cerr << "failed to open " << s_dst_Idx << endl;
	}
	else{
		//*******************************************************
		// first write the beginning of the FileStorage file.
		//*******************************************************
		fs_Idx_write << whatkindImgs << "[:";

		/*each input image*/
		for (int current_image_index = 0; current_image_index != numInputImage; ++current_image_index){
			// save parameters into the ".YML" file for image reconstruction
			// 2-kind Parameters of current image include:
			//     * image size: image_h * image_w;
			//     * possible row and column indices of maximun posterior element of each xPatch in currently read image;
			fs_Idx_write << "{:"
				<< "maxPost_row_col_Idx" << v_maxPost_row_column_Idx[current_image_index]
				<< "}";
		} /* end of each input image*/
	}
	// write the end of the FIleStorage file
	if (!fs_Idx_write.isOpened()){
		cerr << "failed to open " << s_dst_Idx << endl;
	}
	else
		fs_Idx_write << "]";
	fs_Idx_write.release();
}
