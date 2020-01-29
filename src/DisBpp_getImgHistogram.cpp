#include "Dis_Bpp.h"

// /* What does this program do?
// * Loads an image
// * Splits the image into its R, G and B planes using the function split
// * Calculate the Histogram of each 1-channel plane by calling the function calcHist
// * Plot the three histograms in a window
// */

void Distortion_Bpp::getImgHistogram(
	const Mat & src, // only one input image
	const Mat & dst, // output histogram of the only 1 input image
	const int & histSize, // the number of bins , e.g., histSize = 256
	bool & uniform, //
	bool & CurveTypeHist //
	){

	if (!src.data)
	{
		cerr << "No Valid Input Image Being Read Now!\n";
		return;
	}

	/// Separate the image in 3 places ( B, G and R )
	vector<Mat> bgr_planes;
	split(src, bgr_planes);


	/// Set the ranges ( for B,G,R) )
	float range[] = { 0, 256 };
	const float* histRange = { range };

	bool accumulate = false;

	Mat b_hist, g_hist, r_hist;

	// Compute the histograms:
	calcHist(&bgr_planes[0], 1, 0, Mat(), b_hist, 1, &histSize, &histRange, uniform, accumulate);
	calcHist(&bgr_planes[1], 1, 0, Mat(), g_hist, 1, &histSize, &histRange, uniform, accumulate);
	calcHist(&bgr_planes[2], 1, 0, Mat(), r_hist, 1, &histSize, &histRange, uniform, accumulate);

	// Draw the histograms for B, G and R
	int hist_w = 512; int hist_h = 400;
	int bin_w = cvRound((double)hist_w / histSize);

	Mat histImage(hist_h, hist_w, CV_8UC3, Scalar(0, 0, 0));

	/// Normalize the result to [ 0, histImage.rows ]
	normalize(b_hist, b_hist, 0, histImage.rows, NORM_MINMAX, -1, Mat());
	normalize(g_hist, g_hist, 0, histImage.rows, NORM_MINMAX, -1, Mat());
	normalize(r_hist, r_hist, 0, histImage.rows, NORM_MINMAX, -1, Mat());

	// Draw for each channel
	for (int i = 1; i < histSize; i++)
	{
		line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(b_hist.at<float>(i - 1))),
			Point(bin_w*(i), hist_h - cvRound(b_hist.at<float>(i))),
			Scalar(255, 0, 0), 2, 8, 0);
		line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(g_hist.at<float>(i - 1))),
			Point(bin_w*(i), hist_h - cvRound(g_hist.at<float>(i))),
			Scalar(0, 255, 0), 2, 8, 0);
		line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(r_hist.at<float>(i - 1))),
			Point(bin_w*(i), hist_h - cvRound(r_hist.at<float>(i))),
			Scalar(0, 0, 255), 2, 8, 0);
	}

	/// Display
	namedWindow("calcHist Demo", CV_WINDOW_AUTOSIZE);
	imshow("calcHist Demo", histImage);

	waitKey(0);

}

// Please pay attention to the input image -- Mat src.
// It should be of the <uchar> type.
tuple<int, int, int, double, double> Distortion_Bpp::getGrayImgHistogram(
	const Mat & d_src, // only one input image
	const int & histSize, // the number of bins , e.g., histSize = 256 
	const int & hist_w, // e.g., hist_w = 512 or 1024, width of histogram drawing
	const int & hist_h, // e.g., hist_h = 400; height of histogram drawing
	bool  & IsCurveTypeHist, // Flag indicating to draw the histogram, in the form of curves or rectangles
	const string & pathHist, // save histogram
	const string & foutHistVal
	){
#ifdef _DEBUG
	cout << "foutHistVal File Name = " << foutHistVal << "\n";
#endif
	std::ofstream fout_hist;
	fout_hist.open(foutHistVal, std::ofstream::out | std::ofstream::app);

	if (!fout_hist.is_open()){
		std::cout << "File Not Opened" << endl;
	}

	if (!d_src.data){
	
		cerr << "No Valid Input Image Being Read Now!\n";
	}

	// to find the maximum value, usually it is 255
	double maxVal = 0;
	double minVal = 0;
	minMaxLoc(d_src, &minVal, &maxVal, 0, 0);

	// Separate the image in 3 places ( B, G and R )
	// vector<Mat> bgr_planes;
	// split(src, bgr_planes);

	// Compute the histograms:
	bool uniform = true; // Flag indicating whether the histogram is uniform
	bool accumulate = false;
	
	// Set the ranges of intensity
	// e.g., errorImgBit = 8;
	float range[] = { 0, (float)pow(2.0, errorImgBit) };
	const float* histRange = { range };
	Mat grayImg_hist; // it is a matrix with only 1 row but with multiple columns.
	// OpenCV function : calcHist(const Mat* images, ...), where images 每 Source arrays.
	// They all should have the same depth, CV_8U or CV_32F , and the same size. Each of them can have an arbitrary number of channels.
	// So we have to make sure the data-type of the input image -- Mat src -- be of CV_8U or CV_32F type.
	// Changing the data-type of a Mat class instance in OpenCV C++ Interface
	// see http://stackoverflow.com/questions/3188352/changing-the-dataype-of-a-mat-class-instance-in-opencv-c-interface
	Mat src;
	d_src.convertTo(src, CV_32FC1);
	calcHist(&src, 1, 0, Mat(), grayImg_hist, 1, &histSize, &histRange, uniform, accumulate);

	// save the maxVal of the frequency before normalization
	double maxFrequenVal = 0;
	minMaxLoc(grayImg_hist, 0, &maxFrequenVal, 0, 0);
	fout_hist << "Maximum Error Value is " << maxVal << endl;
	fout_hist << "Minimum Error Value is " << minVal << endl;
	fout_hist << "Unnormalized Maximal Frequency is " << maxFrequenVal << endl;


	// get the mean and standard deviation before do the normalization to histogram
	double hist_mean = 0,
		hist_variance = 0;
	int total = std::max<int>(d_src.rows * d_src.cols, 1);
	for (int i = 0; i != histSize; ++i){
		hist_mean += i* grayImg_hist.at<float>(i);
		hist_variance += i * i * grayImg_hist.at<float>(i);
	}

	hist_mean /= total;
	hist_variance /= total;
	hist_variance -= hist_mean * hist_mean;
	cout << "mean = " << hist_mean << ", standard deviation = " << sqrt(hist_variance) << endl;

	// Draw the histograms for gray image
	int bin_w = cvRound((double)hist_w / histSize);

	Mat histImage(hist_h, hist_w, CV_8UC1, Scalar(0));

	/// Normalize the result to [ 0, histImage.rows ]
	// /* function: cv::normalize, to normalize the norm or value range of an array.
	// * C++: void normalize(InputArray src, OutputArray dst, double alpha=1, double beta=0, int norm_type=NORM_L2, int dtype=-1, InputArray mask=noArray() )
	// * Parameters: 
	// *   src 每 input array.
	// *   dst 每 output array of the same size as src .
	// *   alpha 每 norm value to normalize to or the lower range boundary in case of the range normalization.
	// *   beta 每 upper range boundary in case of the range normalization; it is not used for the norm normalization.
	// *   normType 每 normalization type (see the details below).
	// *	dtype 每 when negative, the output array has the same type as src; otherwise, it has the same number of channels as src and the depth = CV_MAT_DEPTH(dtype).
	// *	mask 每 optional operation mask.
	// */

	cv::normalize(grayImg_hist, grayImg_hist, 0, histImage.rows, NORM_MINMAX, -1, Mat());
	cv::Point maxFrequenPoint;
	minMaxLoc(grayImg_hist, 0, &maxFrequenVal, 0, &maxFrequenPoint);
	fout_hist << "Normalized Maximal Frequency is " << maxFrequenVal << endl;
	fout_hist << "Normalized Histogram : Height " << hist_h << " by Width " << hist_w << endl;

	
	// Draw for each channel
	if (IsCurveTypeHist){ // to draw "lines"
		for (int i = 1; i < histSize; i++)
		{
			int tempFrequency = cvRound(grayImg_hist.at<float>(i - 1));

			line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(grayImg_hist.at<float>(i - 1))),
				Point(bin_w*(i), hist_h - cvRound(grayImg_hist.at<float>(i))),
				Scalar(255, 0, 0), 2, 8, 0);
			//		line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(g_hist.at<float>(i - 1))),
			//			Point(bin_w*(i), hist_h - cvRound(g_hist.at<float>(i))),
			//			Scalar(0, 255, 0), 2, 8, 0);
			//		line(histImage, Point(bin_w*(i - 1), hist_h - cvRound(r_hist.at<float>(i - 1))),
			//			Point(bin_w*(i), hist_h - cvRound(r_hist.at<float>(i))),
			//			Scalar(0, 0, 255), 2, 8, 0);
			if (tempFrequency > errHistThres){
				// public member function -- <ios> <iostream> std::ios_base::width
				fout_hist << std::right << "For Bin -- " << i - 1 << " , frequency = " << tempFrequency << endl;
			}
		}
	}

	else{ // to draw "rectangles"
		for (int i = 1; i < histSize; i++){

			int tempFrequency = cvRound(grayImg_hist.at<float>(i - 1));
			rectangle(histImage, Point(bin_w*(i - 1), hist_h),
				Point(bin_w*(i), hist_h - tempFrequency),
				Scalar::all(i - 1),
				CV_FILLED);
			if (tempFrequency > errHistThres){
				fout_hist << std::right <<  "For Bin -- " << i - 1 << " , frequency = " << tempFrequency << endl;
			}

		}
	}
	
	fout_hist.close();

	// get the filtered max and min intensity, i.e., the range of the filtered intensity.
	// filtered means the result after applying the threshold -- errHistThres;
	int minFilterdIntensity = 0;
	int maxFilterdIntensity = 0;

	// /*
	// * break Statement (C++): 
	// * The break statement is used with the conditional switch statement and with the do, for, and while loop statements.
	// * In a switch statement, the break statement causes the program to execute the next statement outside the switch statement. 
	// * Without a break statement, every statement from the matched case label to the end of the switch statement, including the default clause, is executed.
	// * In loops, the break statement ends execution of the nearest enclosing do, for, or while statement. Control passes to the statement that follows the ended statement, if any.
	// * Within nested statements, the break statement ends only the do, for, switch, or while statement that immediately encloses it. 
	// * You can use a return or goto statement to transfer control from more deeply nested structures.
	// */

	for (int i = 0; i != histSize; ++i){
		int tempFrequency = cvRound(grayImg_hist.at<float>(i));
		if (tempFrequency <= errHistThres)
		minFilterdIntensity = i;
		else
			break;
		}

	for (int i = histSize - 1; i > 0; --i){
		int tempFrequency = cvRound(grayImg_hist.at<float>(i));
		if (tempFrequency <= errHistThres)
			maxFilterdIntensity = i;
		else
			break;
	}


	/// Display
	// namedWindow("calcHist Demo", CV_WINDOW_AUTOSIZE);
	// imshow("calcHist Demo", histImage);
#ifdef _DEBUG
	cout << "pathHist = " << pathHist << endl;
#endif
	imwrite(pathHist, histImage);

	tuple<int, int , int, double , double> max_min_filtered_intensity;
	get<0>(max_min_filtered_intensity) = minFilterdIntensity;
	get<1>(max_min_filtered_intensity) = maxFrequenPoint.y; 
	get<2>(max_min_filtered_intensity) = maxFilterdIntensity;
	get<3>(max_min_filtered_intensity) = hist_mean;
	get<4>(max_min_filtered_intensity) = sqrt(hist_variance);
	
	return max_min_filtered_intensity;
}