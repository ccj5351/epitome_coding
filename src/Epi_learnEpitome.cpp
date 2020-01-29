#include "ImageEpitome.hpp"

// this function will just learn one category of images, but with different kinds of parameter setup;
// like to learn the category of Edwin of the Stevens dataset, we collected on the Stevens campus,
// with different kinds of parameters, e.g., epitome size = 128 * 128, 256* 256, or 512 * 512;
void ImgEpitome::learnEpitome(EpitomeVariableName flag){

	////////////////////////////////////
	// Step - Learn Epitome
	////////////////////////////////////

	// input images' names in the current category which will be learned by this function;
	vector<string> filelist;
	GetFileList(DatabaseDir, &filelist);

	// write filelist into file storage, i.e., to save all the images names;
	FileStorage fs_filelist_in(EpitomeResultDir + "/" + whatkindImgs + "_imgNames.yml", FileStorage::WRITE);
	if (!fs_filelist_in.isOpened()){
		cerr << "failed to open " << EpitomeResultDir + "/" + whatkindImgs + "_imgNames.yml" << endl;
		fileStorageHelp();
	}

	fs_filelist_in << whatkindImgs << filelist;
	fs_filelist_in.release();

	// randomly initialize v_eMean
	// vector<double>  v_eMean(eHeight * eWidth, .0), // zeros
	//	v_eVar(eHeight* eWidth, 1.0); // ones
	// Mat temp_img = imread(DatabaseDir + "/" + categories[categoryIdx] + "/" + filelist[0], 0);
	// randomizeEpitomeMean(temp_img, v_eMean);


	switch (flag)
	{
	case ORIGINAL_FFTW: // the 1.0 version of the convolution funciton using FFTW;
		learnImgGenreEpitome(filelist);
		break;
	case UPDATED_FFTW: // the 2.0 version of the convolution funciton using FFTW;
		learnImgGenreEpitome_plus(filelist);
		break;
	default:
		cout << "ERROR! No appropriate parameters was input.\n";
	}
	
	// clear the file name for next I/O operation
	// filelist.clear(); 
	vector<string>().swap(filelist);
}