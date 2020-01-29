#include "filesStorage.hpp"

 void fileStorageHelp()
{
	cout << "\nfilestorage_sample demonstrate the usage of the OpenCV serialization functionality.\n"
		<< "usage:\n"
		<< "output file above can have many different extensions, see below."
		<< "\nThis program demonstrates the use of FileStorage for serialization, that is use << and >>  in OpenCV\n"
		<< "For example, how to create a class and have it serialize, but also how to use it to read and write matrices.\n"
		<< "FileStorage allows you to serialize to various formats specified by the file end type."
		<< "\nYou should try using different file extensions.(e.g. yaml yml xml xml.gz yaml.gz etc...)\n" << endl;
}

 void writeDataToFileStorage(vector<vector<vector<int>>> & v,
	string & filename, // e.g. filename = "test.yml"
	string & fileNodeName // e.g. "imageIdx"
	){
	FileStorage fs_in(filename, FileStorage::WRITE);
	if (!fs_in.isOpened())
	{
		cerr << "failed to open " << filename << endl;
		fileStorageHelp();
	}

	fs_in << fileNodeName << "[";
	for (int i = 0, size1 = v.size(); i < size1; i++){
		// 0: "image height & width for image"
		// 1:  "row indices with random shifting for image"
		// 2: "column indices with random shifting for image" 
		// 3: "row indices of maximum posterior element for image"
		// 4:  "column indices of maximum posterior element for image"
		fs_in << "{:" << "size" << v[i][0] << "row_Idx" << v[i][1]
			<< "col_Idx" << v[i][2] << "maxPost_row_Idx" << v[i][3]
			<< "maxPost_col_Idx" << v[i][4] << "}";
	}
	fs_in << "]";
	fs_in.release();
}



 void readDataFromFileStorage(vector<vector<vector<int>>> & v,
	string & filename, // e.g. filename = "test.yml"
	string & fileNodeName // e.g. "imageIdx"
	){
	FileStorage fs_out(filename, FileStorage::READ);
	if (!fs_out.isOpened())
	{
		cerr << "failed to open " << filename << endl;
		fileStorageHelp();
	}

	for (int i = 0, size1 = v.size(); i < size1; i++){
		FileNode forImages = fs_out[fileNodeName];
		FileNodeIterator it = forImages.begin(), // iterator pointing to the first node element.
			it_end = forImages.end(); // iterator pointing to the element following the last node element.

		// iterate through a sequence using FileNOdeIterator
		for (; it != it_end; it++){
			// 0: "size" = "image height & width for image"
			// 1:  "row_Idx" = "row indices with random shifting for image"
			// 2: "col_Idx" = "column indices with random shifting for image" 
			// 3:  "maxPost_row_Idx" = "row indices of maximum posterior element for image"
			// 4:   "maxPost_col_Idx" = "column indices of maximum posterior element for image"
			(*it)["size"] >> v[i][0];
			(*it)["row_Idx"] >> v[i][1];
			(*it)["col_Idx"] >> v[i][2];
			(*it)["maxPost_row_Idx"] >> v[i][3];
			(*it)["maxPost_col_Idx"] >> v[i][4];
		}
	}
	fs_out.release();
}


void read_Write_FileStorage(vector<vector<vector<int>>> & v,
	string & filename, // e.g. filename = "test.yml"
	string & filenodename,
	int flag  // #define Read_FileStorage 1; #define Write_FileStorage 0
	)
{
	switch (flag){

	case  Write_FileStorage: { //write date to file storage
		writeDataToFileStorage(v, filename, filenodename);
		break;

	}

	case Read_FileStorage: {//read data from file storage
		readDataFromFileStorage(v, filename, filenodename);
		break;
	}

	default:{
		fileStorageHelp();
		cout << "Please input integer 0 for writing or 1 for reading."
			<< "\nother numbers are invalid.\n"; }
	}
}