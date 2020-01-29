#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

/**
* @ see reference http://stackoverflow.com/questions/5089927/show-multiple-2-3-4-images-in-the-same-window-in-opencv
* @brief makeCanvas Makes composite image from the given images
* @param vecMat Vector of Images.
* @param windowHeight The height of the new composite image to be formed.
* @param nRows Number of rows of images. (Number of columns will be calculated
*              depending on the value of total number of images).
* @return new composite image.
*/
cv::Mat makeCanvas(std::vector<cv::Mat>& v_Mat,// Vector of Images
	int windowHeight, // The height of the new composite image to be formed
	int nRows // Number of rows of images. (Number of columns will be calculated
	// depending on the value of total number of images).
	);