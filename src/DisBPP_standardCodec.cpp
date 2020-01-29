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

void Distortion_Bpp::imgCompressViaStandardCodec(
	const string & GetJPEGCompressDir, // = "E:\\ComputerVisionCCJ\\Epitome\\imageDatabase\\simi-objects-bmpImages";
	const string & JPEG2K_EXE_Base_Path
	){

	int imgreadFlag = 0; // 0 = gray, 1 = color images
	// make directories for saving the final compressed reconstruction images
	MakeDir(GetJPEGCompressDir);

	//  check the string is 1-, 2-, or 3-character length, and to make them be 3-character length, if less than 3.
	//    i.e., change '5' to '005', and '15' to '015'; 
	//    Usage: std::string str ("Test string");
	//    std::cout << "The size of str is " << str.length() << " bytes.\n";
	// 
	vector<int> v_jpegParams(2);
	switch (encodeType)
	{
	case JPEG:
		// For JPEG, it can be a quality ( CV_IMWRITE_JPEG_QUALITY ) from 0 to 100 (the higher is the better). Default value is 95.
		v_jpegParams[0] = CV_IMWRITE_JPEG_QUALITY;
		v_jpegParams[1] = jpeg_quality ;
		break;
	case JPEG2000:
		break;
	case PNG:
		// For PNG, it can be the compression level ( CV_IMWRITE_PNG_COMPRESSION ) from 0 to 9. 
		// A higher value means a smaller size and longer compression time. Default value is 3.
		v_jpegParams[0] = CV_IMWRITE_PNG_COMPRESSION;
		v_jpegParams[1] = jpeg_quality % 9;
		break;
	case BMP:
		break;
	default:
		cout << "ERROR! Check the right name of encodeType!\n";
	}

	string s_jpeg_q = std::to_string(static_cast<long long>(jpeg_quality));
	if (s_jpeg_q.length() == 1){
		s_jpeg_q = "00" + s_jpeg_q;
	}

	if (s_jpeg_q.length() == 2){
		s_jpeg_q = "0" + s_jpeg_q;
	}


	// to get the file lists of each input image category,
	vector<string> filelist;
	// get image names in each category
	    GetFileList(DatabaseDir, &filelist);
		int fileListSize = filelist.size();
		std::string nameTemp = GetJPEGCompressDir + "/" + "cq-" + s_jpeg_q;
		MakeDir(nameTemp);

		for (int fileIdx = 0; fileIdx != fileListSize; ++fileIdx){
			// read input image
			cv::Mat input = imread(DatabaseDir + "/" + filelist[fileIdx], imgreadFlag); // read gray or color images
			
			// std::size_t pos = filelist[fileIdx].find(".");
			// string imageNameNoFileExtension = filelist[fileIdx].substr(0, pos);
			string nameTemp_J2k = nameTemp + "/" + "J2K";
			if (encodeType == JPEG2000)
				MakeDir(nameTemp_J2k);

			string img_save_path = nameTemp + "/" +  filelist[fileIdx];
			setFileExtension(img_save_path, encodeTypeforErrorImgs);
			
			if (encodeType == BMP){
				imwrite(img_save_path, input);
			}

			else if (encodeType == JPEG2000){
				char * p_decompress_buffer = new char[MAX_CHAR_NUM_OF_FILES_PATH];
				char * p_compress_buffer = new char[MAX_CHAR_NUM_OF_FILES_PATH];

				// syntax: char * strcpy(char * destination, const char * source);
				// Copies the C string pointed by source into the array pointed by destination, 
				// including the terminating null character (and stopping at that point).
				string img_path = DatabaseDir + "/" + filelist[fileIdx];
				char * char_img_path = new char[img_path.length() + 1];
				std::strcpy(char_img_path, img_path.c_str());

				string de_img_path = nameTemp + "/" + nameDifference + filelist[fileIdx];;
				char * char_de_img_save_path = new char[de_img_path.length() + 1];
				std::strcpy(char_de_img_save_path, de_img_path.c_str());

				string j2k_img_save_path = nameTemp_J2k + "/" + filelist[fileIdx];
				setFileExtension(j2k_img_save_path, encodeTypeforErrorImgs);
				char * char_j2k_img_save_path = new char[j2k_img_save_path.length() + 1];
				std::strcpy(char_j2k_img_save_path, j2k_img_save_path.c_str());


				string JPEG2K_Compress_Path = JPEG2K_EXE_Base_Path + "/" + "opj_compress.exe";
				string JPEG2K_Decompress_Path = JPEG2K_EXE_Base_Path + "/" + "opj_decompress.exe";
				int J2K_Compression_Ratio = 101 - jpeg_quality; // jpeg_quality belongs to [0, 100]
				// that means J2K_Compression_Ratio belongs to [1, 101];
				string s_J2K_Compression_Ratio = to_string(static_cast<long long>(J2K_Compression_Ratio));

				char * char_J2K_Compression_Ratio = new char[s_J2K_Compression_Ratio.length() + 1];
				std::strcpy(char_J2K_Compression_Ratio, s_J2K_Compression_Ratio.c_str());
				char * char_JPEG2K_Compress_Path = new char[JPEG2K_Compress_Path.length() + 1];
				std::strcpy(char_JPEG2K_Compress_Path, JPEG2K_Compress_Path.c_str());
				char * char_JPEG2K_Decompress_Path = new char[JPEG2K_Decompress_Path.length() + 1];
				std::strcpy(char_JPEG2K_Decompress_Path, JPEG2K_Decompress_Path.c_str());

				sprintf(p_compress_buffer, "%s -i %s -o %s -r %s", char_JPEG2K_Compress_Path, char_img_path, char_j2k_img_save_path, char_J2K_Compression_Ratio);
				sprintf(p_decompress_buffer, "%s -i %s -o %s", char_JPEG2K_Decompress_Path, char_j2k_img_save_path, char_de_img_save_path);

				// compress, from .bmp to .j2k;
				std::system(p_compress_buffer);
				// decompress, from .jk2 to .bmp;
				std::system(p_decompress_buffer);
				
				delete[] char_img_path;
				delete[] char_de_img_save_path;
				delete[] char_j2k_img_save_path;
				delete[] char_J2K_Compression_Ratio;
				delete[] char_JPEG2K_Compress_Path;
				delete[] char_JPEG2K_Decompress_Path;
				delete[] p_decompress_buffer;
				delete[] p_compress_buffer;
				cerr << "Currently JPEG2000 is not supported! Please try other methods!\n";
			}
			else
				imwrite(img_save_path, input, v_jpegParams);
		} /*end of each image*/
		
		// release memory
		vector<string>().swap(filelist);
}

