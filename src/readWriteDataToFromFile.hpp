#ifndef _read_Write_Data_To_From_File_HPP
#define _read_Write_Data_To_From_File_HPP
#include <iostream>
#include <fstream>
#include <string>
#include <opencv2/core/core.hpp>
using namespace std;
using namespace cv;

void readDoubleDataFromBinaryFile(string filename, // "E:\\ComputerVisionCCJ\\Epitome\\eMean.bin"
	std::vector<double> & v_dataFromFile);

void readDoubleDataFromTxtFile(string filename, // e.g. "eMean.txt"
	std::vector<double> & v_dataFromFile);


void writedoubleMatToTxtFile(cv::Mat_<double> & m,
	string filename // e.g. "eMean.txt"
	);


// save the index to txt files
void writeIdxVectorDataToTxtFile(vector<vector<vector<int>>> & v, string filename);
#endif