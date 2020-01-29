#include "readWriteDataToFromFile.hpp"

void readDoubleDataFromBinaryFile(string filename, // "E:\\ComputerVisionCCJ\\Epitome\\eMean.bin"
	std::vector<double> & v_dataFromFile){
	long length;
	ifstream readFile;
	//	string filename = "E:\\ComputerVisionCCJ\\Epitome\\eMean.bin";
	readFile.open(filename, ios::binary | ios::in);

	// get length of file:
	// "seekg" is used to move the position to the end of the file, and then back to the beginning.
	readFile.seekg(0, ios::end);
	length = readFile.tellg();
	readFile.seekg(0, ios::beg);
	if (!readFile){
		cerr << "Open error!" << endl;
		//		exit(1);
	}

	if (abs(length) > (v_dataFromFile.size() * sizeof(double))){
		cout << "sizeof(double) = " << sizeof(double) << endl;
		cerr << "The vector is not big enough to save all the data from the binary file.\n";
		//		exit(1);
	}

	// read data as a block:
	readFile.read(reinterpret_cast<char*>(&v_dataFromFile[0]), length);
	readFile.close();
	cout << " The binary file is of length " << length << " bytes, i.e., " << length / sizeof(double) << " doubles\n";
	cout << "The output vector data: its size is " << v_dataFromFile.size() << endl;
	cout << "The first two elements: " << v_dataFromFile[0] << " and  " << v_dataFromFile[1] << endl;
	cout << "The last two elements: " << v_dataFromFile[length / sizeof(double) - 2] << " and " << v_dataFromFile[length / sizeof(double) - 1] << endl;
}

void readDoubleDataFromTxtFile(string filename, // e.g. "eMean.txt"
	std::vector<double> & v_dataFromFile){
	// open file    
	ifstream inputFile(filename);

	// test file open   
	if (inputFile) {
		double value;
		// read the elements in the file into a vector  
		while (inputFile >> value) {
			v_dataFromFile.push_back(value);
		}
	}

	inputFile.close();// close the file
}


void writedoubleMatToTxtFile(cv::Mat_<double> & m,
	string filename // e.g. "eMean.txt"
	){

	ofstream fout(filename);

	if (!fout){
		cout << "File Not Opened" << endl;  return;
	}

	for (int i = 0, rows = m.rows; i < rows; i++){
		for (int j = 0, cols = m.cols; j < cols; j++){
			fout << m.at<double>(i, j) << " ";
		}
		fout << endl;
	}
	fout.close();
}


// save the index to txt files
void writeIdxVectorDataToTxtFile(vector<vector<vector<int>>> & v, string filename){

	ofstream fout(filename);
	if (!fout){
		std::cout << "File Not Opened" << endl;  return;
	}

	for (int i = 0, size1 = v.size(); i < size1; i++){

		for (int j = 0; j < 5; j++){
			switch (j){
			case 0: {
				fout << "image height & width for image" << i + 1 << "/" << size1 << ":\n";
				for (int k = 0, length = v[i][j].size(); k < length; k++){
					fout << v[i][j][k] << " ";
				}
				fout << endl;
				break;
			}
			case 1: {
				fout << "row indices with random shifting for image" << i + 1 << "/" << size1 << ":\n";
				for (int k = 0, length = v[i][j].size(); k < length; k++){
					fout << v[i][j][k] << " ";
				}
				fout << endl;
				break;
			}
			case 2: {
				fout << "column indices with random shifting for image" << i + 1 << "/" << size1 << ":\n";
				for (int k = 0, length = v[i][j].size(); k < length; k++){
					fout << v[i][j][k] << " ";
				}
				fout << endl;
				break;
			}
			case 3: {
				fout << "row indices of maximum posterior element for image" << i + 1 << "/" << size1 << ":\n";
				for (int k = 0, length = v[i][j].size(); k < length; k++){
					fout << v[i][j][k] << " ";
				}
				fout << endl;
				break;
			}

			default: {
				fout << "column indices of maximum posterior element for image" << i + 1 << "/" << size1 << ":\n";
				for (int k = 0, length = v[i][j].size(); k < length; k++){
					fout << v[i][j][k] << " ";
				}
				fout << endl;
			}
			}
		} // end of j
		fout << endl;
		fout << endl;
	}// end of i

	fout.close();
}