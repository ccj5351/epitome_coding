#include "ImageEpitome.hpp"

#define _CRT_SECURE_NO_DEPRECATE
#include <stdio.h>
#include <math.h>
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "iostream"
#include "fstream"
using namespace std;
using namespace cv;
// encode the max-row-column indices
// DPCM -- Differential Pulse Code Modulation on DC components in JPEG
// DC component is large and varied, but often closed to previous value (like lossless JPEG)
// encode the difference from previous value


void ImgEpitome::encodeMaxRowColIdx(){
	// do nothing now
	// to be continued ...
}


// decode the max-row-column indices
void ImgEpitome::decodeMaxRowColIdx(){

	// do nothing now
	// to be continued ...
}

