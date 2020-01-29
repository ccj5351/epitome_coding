#include "ShowManyImages.h"
using namespace std;
using namespace cv;
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
	){
	int N = v_Mat.size();
	nRows = nRows > N ? N : nRows;
	int edgeThickness = 10;
	int imagesPerRow = ceil(double(N) / nRows);
	int resizeHeight = floor(2.0 * ((floor(double(windowHeight - edgeThickness) / nRows)) / 2.0)) - edgeThickness;
	int maxRowLength = 0;

	std::vector<int> resizeWidth;
	for (int i = 0; i < N;){
		int thisRowLen = 0;
		for (int k = 0; k < imagesPerRow; k++) {
			double aspectRatio = double(v_Mat[i].cols) / v_Mat[i].rows;
			int temp = int(ceil(resizeHeight * aspectRatio));
			resizeWidth.push_back(temp);
			thisRowLen += temp;
			if (++i == N) break;
		}
		if ((thisRowLen + edgeThickness * (imagesPerRow + 1)) > maxRowLength){
			maxRowLength = thisRowLen + edgeThickness * (imagesPerRow + 1);
		}
	}
	int windowWidth = maxRowLength;
	cv::Mat canvasImage(windowHeight, windowWidth, CV_8UC1, Scalar(0));

	for (int k = 0, i = 0; i < nRows; i++){
		int y = i * resizeHeight + (i + 1) * edgeThickness;
		int x_end = edgeThickness;
		for (int j = 0; j < imagesPerRow && k < N; k++, j++){
			int x = x_end;
			cv::Rect roi(x, y, resizeWidth[k], resizeHeight);
			cv::Mat target_ROI = canvasImage(roi);
			cv::resize(v_Mat[k], target_ROI, target_ROI.size());
			x_end += resizeWidth[k] + edgeThickness;
		}
	}
	return canvasImage;
}