#include "Dis_Bpp.h"
/*

string DatabaseDir ;
string EpitomeResultDir;
string ReconsCompresDir;
EncodeTypeName encodeType; // e.g., jpeg, jpeg2000, bmp, png, etc.
QuantizeType quantizeType;
string s_quantizeType;
const int errorImgBit ;
const int jpeg_quality;
const string encodeTypeforReconsImgs; // e.g. = ".png";
const string encodeTypeforErrorImgs; // e.g. = ".jpg";
const string nameDifference; // e.g., = "recon-";
const string whatKindofImgs;
const int errHistThres;

*/

//*********************************************************************************************************
// to get rate-distortion data for the JPEG-compressed-images of some compression quality, like 90.
//*********************************************************************************************************
tuple<double, double, double> Distortion_Bpp::getStandardCodecRate_Distortion(
	const string & GetJPEGCompressDir // some compression quality images' directory
	){

	int imgReadFlag = 0; // 0 = gray, 1 = color images
	vector<string> input_filelist;
	GetFileList(DatabaseDir, &input_filelist);
	int fileListSize = input_filelist.size();
	vector<tuple<uintmax_t, uintmax_t, double, double>> v_data(fileListSize, tuple<uintmax_t, uintmax_t, double, double>());
	//	tuple<0> -- input_fileSize;
	//	tuple<1> -- compressed-error-images_fileSize;
	//	tuple<2> -- psnr;
	//	tuple<3> -- bpp;
	uintmax_t sum_input_size = 0;
	uintmax_t sum_errorimg_size = 0;
	string recon_path;

	for (int fileIdx = 0; fileIdx != fileListSize; ++fileIdx){
		// input image's file size (Units : Bytes)
		string input_path = DatabaseDir + "/" + input_filelist[fileIdx];
		get<0>(v_data[fileIdx]) = file_size(input_path);
		sum_input_size += get<0>(v_data[fileIdx]);

		// make sure the filename is imencodeType, say, ".jpg"
		std::size_t pos = input_filelist[fileIdx].find(".");      // e.g., position of ".jpg"
		string imageNameNoFileExtension = input_filelist[fileIdx].substr(0, pos);

		
		string recon_path_j2k = GetJPEGCompressDir + "/" + nameDifference + imageNameNoFileExtension + ".bmp";
		

		if (encodeType == JPEG2000){
			recon_path = GetJPEGCompressDir + "/" + "J2K";
		}
		else{
		recon_path = GetJPEGCompressDir + "/" + nameDifference + input_filelist[fileIdx];
		setFileExtension(recon_path, encodeTypeforErrorImgs);
		// get jpeg compressed image's file size (units : Bytes)
		get<1>(v_data[fileIdx]) = file_size(recon_path);
		sum_errorimg_size += get<1>(v_data[fileIdx]);
		}

		 

		// to get PSNR
		Mat_<double> input = imread(input_path, imgReadFlag);
		Mat_<double> recon;
		if (encodeType == JPEG2000){
		 recon = imread(recon_path_j2k, imgReadFlag);
		}
		else {
			recon = imread(recon_path, imgReadFlag);
		}
		get<2>(v_data[fileIdx]) = getPSNR(input, recon);


		// to get bits per pixel
		int numPixel = input.rows * input.cols;
		get<3>(v_data[fileIdx]) = (double)numPixel;
	} /*end of each image*/

	// to get weighted average psnr 
	// to get weighted average bpp
	if (encodeType == JPEG2000)
		sum_errorimg_size = getFileSizeFromDir(recon_path);
	double weight = 0;
	double average_psnr = 0; // weighted average psnr (peak signal-to-noise ratio)
	double average_bpp = 0; // weighted average bpp (bits per pixel)
	double compressionRatio = (double(sum_input_size)) / (double(sum_errorimg_size)); // original files' size divided by the compressed files' size
	for (int fileIdx = 0; fileIdx != fileListSize; ++fileIdx){
		if (encodeType == JPEG2000)
			weight = 1.0 / fileListSize;
		else
			weight = (double(get<1>(v_data[fileIdx]))) / sum_errorimg_size;
		average_psnr += get<2>(v_data[fileIdx]) * weight;
		average_bpp += get<3>(v_data[fileIdx]);
	}
	average_bpp = 8 * sum_errorimg_size / average_bpp;// bit per pixel
	return tuple<double, double, double>(average_psnr, average_bpp, compressionRatio);
}


void Distortion_Bpp::getStandardCodecRate_DistortionfromDir(
	const string & GetJPEGCompressDir // = "E:\\imageDatabase\\forPSNR\\compress-bmp-images\\compress-small-bmp-images";
	){

	int imgReadFlag = 0; // 0 = gray, 1 = color images
	vector<string> compressQualityCategories, filelist;
	GetDirList(GetJPEGCompressDir, &compressQualityCategories);
	int compressQualityCategorySize = compressQualityCategories.size();
	// for saving the tuples of (psnr, bpp, compression ratio)
	vector<tuple<double, double, double>> v_datas(compressQualityCategorySize);
	string fileStreamName = GetJPEGCompressDir + "/" + whatKindofImgs + "-psnr-bpp-cr.txt";
	ofstream fout_psnr(fileStreamName, std::ofstream::out | std::ofstream::app);

	if (!fout_psnr){
		std::cout << "File Not Opened" << endl;
	}

	for (int CQIdx = 0; CQIdx != compressQualityCategorySize; ++CQIdx){ // for each compression quality
		fout_psnr << whatKindofImgs << "'s " <<compressQualityCategories[CQIdx] + " begin here ----\n";
		string tempPath = GetJPEGCompressDir + "/" + compressQualityCategories[CQIdx];
		tuple<double, double, double> t_output = getStandardCodecRate_Distortion(tempPath);
			fout_psnr <<  " -- weighted average PSNR = " << get<0>(t_output) << endl;
			fout_psnr <<  " -- weighted average bpp = " << get<1>(t_output) << endl;
			fout_psnr <<  " -- Compression Ratio = " << get<2>(t_output) << endl;
			// save the tuple(PSNR, bpp, CR)
			v_datas[CQIdx] = t_output;
		fout_psnr << endl;
	} // // end of each compression quality
	fout_psnr.close();


	//  the txt file, which can be used to send to GNUPlot for plotting
	// to the tuple(PSNR, bpp, CR)
		ofstream fout_psnr_txt( GetJPEGCompressDir + "/"  + whatKindofImgs + "-" + s_encodeType + "-rd.txt");
		if (!fout_psnr_txt){
			std::cout << "File Not Opened" << endl;
		}
		for (int CQIdx = 0; CQIdx != compressQualityCategorySize; ++CQIdx){ // for each compression quality
			fout_psnr_txt << get<0>(v_datas[CQIdx]) << " "
				<< get<1>(v_datas[CQIdx]) << " "
				<< get<2>(v_datas[CQIdx]) << endl;
		}  // end of each compression quality
		fout_psnr_txt.close();
}