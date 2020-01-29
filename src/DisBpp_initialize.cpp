#include "Dis_Bpp.h"

/*
string DatabaseDir ;
string EpitomeResultDir;
string ReconsCompresDir;
EncodeTypeName encodeType; // e.g., jpeg, jpeg2000, bmp, png, etc.
const int errorImgBit ;
const int jpeg_quality;
const string encodeTypeforReconsImgs; // e.g. = ".png";
const string encodeTypeforErrorImgs; // e.g. = ".jpg";
const string nameDifference; // e.g., = "recon-";
*/

// constructor
Distortion_Bpp::Distortion_Bpp( string & input, string & epi, string & compre, string & s_standard_encodeType, 
	EncodeTypeName & encodeType, QuantizeType & quantizetype,
	string & s_quantizeType, int  & errImgBit, int & jpegQual,
	string & encodeTypeforReconsImgs,
	string & encodeTypeforErrorImgs,
	string & nameDifference,
	string & whatKindofImgs,
	int & errHistThres,
	int & save_img_num,
	bool & isGary
	) :
	errorImgBit (errImgBit),
	jpeg_quality(jpegQual),
	encodeTypeforReconsImgs(encodeTypeforReconsImgs),
	encodeTypeforErrorImgs(encodeTypeforErrorImgs),
	nameDifference(nameDifference),
	whatKindofImgs(whatKindofImgs), 
	s_encodeType(s_standard_encodeType),
	errHistThres(errHistThres),
	save_img_number(save_img_num),
	isGrayScale(isGary)
{
	this->DatabaseDir = input;
	this->EpitomeResultDir = epi;
	this->ReconsCompresDir = compre;
	this ->encodeType = encodeType;
	this->quantizeType = quantizetype;
	this->s_quantizeType = s_quantizeType;
}


// get PSNR of 1 image;
double Distortion_Bpp::getPSNR(const Mat& input, const Mat& recons){

	Mat s1;
	absdiff(input, recons, s1);       // |I1 - I2|
	s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
	s1 = s1.mul(s1);           // |I1 - I2|^2

	Scalar s = sum(s1);         // sum elements per channel

	double sse = s.val[0] + s.val[1] + s.val[2]; // sum channels

	if (sse <= 1e-10) // for small values return zero
		return 0;
	else
	{
		double  mse = sse / (double)(input.channels() * input.total());
		double psnr = 10.0*log10((255 * 255) / mse);
		return psnr;
	}
}


void Distortion_Bpp::parseString(const string & s_input, // e.g., == "4-8-12"
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
		std::cout << *(i_output + j) << " ";
	}
	std::cout << endl;
#endif // _DEBUG 


}


