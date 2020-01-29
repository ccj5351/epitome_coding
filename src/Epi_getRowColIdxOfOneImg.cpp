#include "ImageEpitome.hpp"

tuple<vector<int>, vector<int>> ImgEpitome::getRowColIdxOfOneImg(
	const int & image_h, // image height
	const int & image_w // image width
	// const int & patchSideLengh, // the side length of the patch
	// const int & denominator, // e.g., denominator = 2, 
	// bool displayProcess // for display result
	){

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// generating possible patches extracted out of the current input image
	// last possible location of patches in the input image x;

	int c_xEndIdx = image_w - patchSideLengh;
	int r_xEndIdx = image_h - patchSideLengh;
	
	// e.g., denominator = 2, 
	// int patchSpacing = patchSideLengh / denominator;
	// patchSpacing = patchSpacing > 0 ? patchSpacing : 0;

	// possible indices of each xPatch in currently read image
	std::vector<int> v_rIdx(0);
	std::vector<int> v_cIdx(0);

	if (verbose){
		std::cout << "\nc_xEndIdx = " << c_xEndIdx << ", r_xEndIdx = " << r_xEndIdx << "\n";
	}

	int patchWiggle = (patchSpacing / 2) > 1 ? (patchSpacing / 2) : 1 ;
	//	if (patchWiggle == 0) patchWiggle = 1; // avoid mistake for the case of patchSpacing < 2.


	for(int r_Idx = 0; r_Idx < r_xEndIdx; r_Idx += patchSpacing){
		v_rIdx.push_back(r_Idx);	
	}

	for (int c_Idx = 0; c_Idx < c_xEndIdx; c_Idx += patchSpacing){
		v_cIdx.push_back(c_Idx);
	}

	// include the last possible location of patches in the input image x;
	v_rIdx.push_back(r_xEndIdx);
	v_cIdx.push_back(c_xEndIdx);

	if (verbose){
		int r_size = v_rIdx.size(), c_size = v_cIdx.size();
		int PatchNum = r_size* c_size;
		std::cout << "\n The size of v_rIdx is: " << r_size
			<< "\n The size of v_cIdx is: " << c_size
			<< "\n The number of possible patches is: " << PatchNum
			<< endl;
	}

	tuple<vector<int>, vector<int>> tuple_vector_out;
	get<0>(tuple_vector_out) = v_rIdx;
	get<1>(tuple_vector_out) = v_cIdx;

	return tuple_vector_out;

}