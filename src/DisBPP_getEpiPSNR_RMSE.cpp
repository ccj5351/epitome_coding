#include "Dis_Bpp.h"

// In the absence of noise, the two images I and I-recon are identical, 
// and thus the MSE is zero. In this case the PSNR is infinite (or undefined, see Division by zero)
void Distortion_Bpp::getEpiPSNR_MSE(const Mat& input,
	const Mat& recons, double & psnr, double & mse){
	Mat s1;
	absdiff(input, recons, s1);       // |I1 - I2|
	s1.convertTo(s1, CV_32F);  // cannot make a square on 8 bits
	s1 = s1.mul(s1);           // |I1 - I2|^2

	Scalar s = sum(s1);         // sum elements per channel

	double sse = s.val[0] + s.val[1] + s.val[2]; // sum channels

	if (sse <= 1e-10){ // for small values return zero
		psnr = 0;
		mse = 0;
	}
	else
	{
		mse = sse / (double)(input.channels() * input.total());
		psnr = 10.0*log10((255 * 255) / mse);
	}
}


void Distortion_Bpp::getEpiPSNR_MSE_FromDir(string & inputDir,
	string & epi_res_Dir,
	string nameDifference // = "recon-"
	){
	std::cout << "Calculating ... \n";
	string f_name = epi_res_Dir + "/" + "psnr_rmse.txt",
		f_name_2 = epi_res_Dir + "/" + "gnuplot.txt";
	// fstream f_out(f_name, std::fstream::out | std::fstream::app)
	fstream f_out(f_name, std::fstream::out),
		f_out_2(f_name_2, std::fstream::out);

	vector<string> v_imageFileNameList;
	GetFileList(inputDir, &v_imageFileNameList);
	int imgNum = v_imageFileNameList.size();

	vector<double> v_psnr(imgNum, 0.0), v_mse(imgNum, 0.0);

	vector<string> Dirs_epiLearnParam, Dirs_EpiReconParam;
	GetDirList(epi_res_Dir, &Dirs_epiLearnParam);

	int Dirs_epiLearnParam_Num = Dirs_epiLearnParam.size();
	for (int i = 0; i != Dirs_epiLearnParam_Num; ++i){
		string tempDir = epi_res_Dir + "/" + Dirs_epiLearnParam[i];
		Dirs_EpiReconParam.clear();
		GetDirList(tempDir, &Dirs_EpiReconParam);
		int Dirs_EpiReconParam_Num = Dirs_EpiReconParam.size();
		for (int j = 0; j != Dirs_EpiReconParam_Num; ++j){
			for (int imgIdx = 0; imgIdx != imgNum; ++imgIdx){
				Mat inputImg = imread(inputDir + "/" + v_imageFileNameList[imgIdx], img_read_flag);
				Mat reconImg = imread(tempDir + "/" + Dirs_EpiReconParam[j] + "/" + nameDifference + v_imageFileNameList[imgIdx], img_read_flag);
				getEpiPSNR_MSE(inputImg, reconImg, v_psnr[imgIdx], v_mse[imgIdx]);
			} // end of images
			// calculate average psnr and mse

			// initialized value is important here, 
			// here we use a larger number for min_psnr and min_rmse, like "10000";
			// but we use a small number for max_psnr and max_rmse, like "0";
			double average_psnr = 0, average_mse = 0, max_psnr = 0, max_mse = 0,
				   min_psnr = 10000, min_mse = 10000;

			for (std::vector<double>::iterator k = v_psnr.begin(); k != v_psnr.end(); ++k){
				average_psnr += *k;
				max_psnr = (max_psnr >= *k) ? max_psnr : *k;
				min_psnr = (min_psnr <= *k) ? min_psnr : *k;
			}
			average_psnr /= imgNum;

			for (std::vector<double>::iterator k = v_mse.begin(); k != v_mse.end(); ++k){
				average_mse += *k;
				max_mse = (max_mse >= *k) ? max_mse : *k;
				min_mse = (min_mse <= *k) ? min_mse : *k;
			}
			average_mse /= imgNum;

			// for testing 
			auto min_max_psnr_idx = std::minmax_element(v_psnr.begin(), v_psnr.end());
			f_out << "min psnr at: " << (min_max_psnr_idx.first - v_psnr.begin())
				<< ", = " << *min_max_psnr_idx.first
				<< ", i.e., image " << v_imageFileNameList[min_max_psnr_idx.first - v_psnr.begin()]
				<< '\n';

			f_out << "max psnr at: " << (min_max_psnr_idx.second - v_psnr.begin())
				<< ", = " << *min_max_psnr_idx.second
				<< ", i.e., image " << v_imageFileNameList[min_max_psnr_idx.second - v_psnr.begin()]
				<< '\n';

			auto min_max_mse_idx = std::minmax_element(v_mse.begin(), v_mse.end());
			f_out << "min mse at: " << (min_max_mse_idx.first - v_mse.begin())
				<< ", = " << *min_max_mse_idx.first
				<< ", i.e., image " << v_imageFileNameList[min_max_mse_idx.first - v_mse.begin()]
				<< '\n';
			f_out << "max mse at: " << (min_max_mse_idx.second - v_mse.begin())
				<< ", = " << *min_max_mse_idx.second
				<< ", i.e., image " << v_imageFileNameList[min_max_mse_idx.second - v_mse.begin()]
				<< '\n';


			f_out << Dirs_epiLearnParam[i] << "/" << Dirs_EpiReconParam[j] << endl
				<< "average psnr = " << average_psnr << ", average mse =" << average_mse
				<< ", min psnr = " << min_psnr << ", max psnr = " << max_psnr
				<< ", min mse = " << min_mse << ", max mse = " << max_mse
				<< endl << endl;
			f_out_2 << average_psnr << " " << average_mse << " ";

		} // end of Dirs_EpiReconParam loop
		f_out_2 << endl;
	} // end of Dirs_epiLearnParam loop
	// release memory
	vector<double>().swap(v_psnr);
	vector<double>().swap(v_mse);
	vector<string>().swap(v_imageFileNameList);
	vector<string>().swap(Dirs_epiLearnParam);
	vector<string>().swap(Dirs_EpiReconParam);
	std::cout << "Function getPSNR_MSE_FromDir Finishes!\n";
}