#include "Dis_Bpp.h"
/*
**************************************************************************
*** Bit - plane slicing
**************************************************************************
*/

cv::Mat Distortion_Bpp::get_i_th_BitPlane(const e_bit_plane & e_bp, s_bit_plane & s_bp, const cv::Mat & in_img){
	switch (e_bp){
	case plane_8:
		s_bp.i_th_bit = 8; s_bp.s_bit_plane = "plane_8";
		break;
	case plane_7:
		s_bp.i_th_bit = 7; s_bp.s_bit_plane = "plane_7";
		break;
	case plane_6:
		s_bp.i_th_bit = 6; s_bp.s_bit_plane = "plane_6";
		break;
	case plane_5:
		s_bp.i_th_bit = 5; s_bp.s_bit_plane = "plane_5";
		break;
	case plane_4:
		s_bp.i_th_bit = 4; s_bp.s_bit_plane = "plane_4";
		break;
	case plane_3:
		s_bp.i_th_bit = 3; s_bp.s_bit_plane = "plane_3";
		break;
	case plane_2:
		s_bp.i_th_bit = 2; s_bp.s_bit_plane = "plane_2";
		break;
	case plane_1:
		s_bp.i_th_bit = 1; s_bp.s_bit_plane = "plane_1";
		break;
	default:
		cout << "ERROR! NO appropriate type is input!\n";
	}

	s_bp.mask = 1 << (s_bp.i_th_bit - 1);
	cout << s_bp.s_bit_plane << endl;
	cv::Mat  out_img = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
	// cout << out_img.at<uchar>(0, 0);
	// cout << in_img.at<uchar>(0, 0);
	for (int i = 0, nrows = in_img.rows, ncols = in_img.cols; i < nrows; i++)
		for (int j = 0; j < ncols; j++){    
			out_img.at<uchar>(i, j) = ((in_img.at<uchar>(i, j) & s_bp.mask) == 0 ? 0 : 1);
		}
	return out_img;
}

cv::Mat Distortion_Bpp::recon_via_4_MSB_BitPlane(const cv::Mat & in_img){ // reconstruct images via the 4 most significant important bit-planes
	e_bit_plane e_bp_8 = plane_8, e_bp_7 = plane_7;
	e_bit_plane e_bp_6 = plane_6, e_bp_5 = plane_5;
	/*e_bit_plane e_bp_4 = plane_4, e_bp_3 = plane_3;
	e_bit_plane e_bp_2 = plane_2, e_bp_1 = plane_1;
	*/
	s_bit_plane s_bp_8, s_bp_7, s_bp_6, s_bp_5;
	/*s_bit_plane s_bp_4, s_bp_3, s_bp_2, s_bp_1;*/

	cv::Mat out_img_8, out_img_7, out_img_6, out_img_5;
	/*cv::Mat  out_img_4, out_img_3, out_img_2, out_img_1, out_img;
	*/
	if (!in_img.empty()){
		out_img_8 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
		out_img_7 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
		out_img_6 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
		out_img_5 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
		//	out_img_4 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
		//	out_img_3 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
		//	out_img_2 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
		//	out_img_1 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
	}
	else
		std::cout << "ERROR! No image has been read!\n";

	out_img_8 = get_i_th_BitPlane(e_bp_8, s_bp_8, in_img);
	out_img_7 = get_i_th_BitPlane(e_bp_7, s_bp_7, in_img);
	out_img_6 = get_i_th_BitPlane(e_bp_6, s_bp_6, in_img);
	out_img_5 = get_i_th_BitPlane(e_bp_5, s_bp_5, in_img);
	//	get_i_th_BitPlane(e_bp_4, s_bp_4, in_img, out_img_4);
	//	get_i_th_BitPlane(e_bp_3, s_bp_3, in_img, out_img_3);
	//	get_i_th_BitPlane(e_bp_2, s_bp_2, in_img, out_img_2);
	//	get_i_th_BitPlane(e_bp_1, s_bp_1, in_img, out_img_1);
	cv::Mat  out_img = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
	out_img += out_img_5 * 16 + out_img_6 * 32 + out_img_7 * 64 + out_img_8 * 128;
	/*out_img = out_img_1 * 255 * coef[0] + out_img_2 * 255 * coef[1] + out_img_3 * 255 * coef[2] + out_img_4 * 255 * coef[3] +
	out_img_5 * 255 * coef[4] + out_img_6 * 255 * coef[5] + out_img_7 * 255 * coef[6] + out_img_8 * 255 * coef[7];*/
	
	return out_img;

}

cv::Mat Distortion_Bpp::recon_via_2_MSB_BitPlane(const cv::Mat & in_img){ // reconstruct images via the 2 most significant important bit-planes
	e_bit_plane e_bp_8 = plane_8, e_bp_7 = plane_7;
	s_bit_plane s_bp_8, s_bp_7;
	cv::Mat out_img_8, out_img_7;
	if (!in_img.empty()){
		out_img_8 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
		out_img_7 = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
	}
	else
		std::cout << "ERROR! No image has been read!\n";

	out_img_8 = get_i_th_BitPlane(e_bp_8, s_bp_8, in_img);
	out_img_7 = get_i_th_BitPlane(e_bp_7, s_bp_7, in_img);
	cv::Mat  out_img = cv::Mat(in_img.size(), CV_8UC1, Scalar(0));
	out_img +=  out_img_7 * 64 + out_img_8 * 128;
	return out_img;
}
