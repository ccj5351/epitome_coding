#ifndef _read_Write_FileStorage_HPP
#define _read_Write_FileStorage_HPP
// you can store and then restore various Opencv data structures ro/from XML or YAML formats.
// Also, it is possible to store and load arbitrary complex data structures, which include OpenCV data structures,
// as well as primitive data types (integer ans floating-point numbers and text strings) as their elements.

//using the following procedure to write something to XML or YAML
// 1. create new FileStorage and open it for writing.  It can be done with a single call to FileStorage::FileStorage() constructor,
// that takes a filename, or you can use the default constructor and then call FileStorage::open(). Format of the file (XML or YAML) is determined from
// the filename extension (".xml" and ".yml"/".yaml", respectively)
// 2. Write all the data you want using the streaming operator << just like in the case of STL streams.
// 3. Close the file using FileStorage::release(). FileStorage destructor also closes the file.

// using the following procedure to read something from XML or YAML
// 1. Open the file storage using FileStorage::FileStorage() constructor or FileStorage::open() method.
//    In the current implementation the whole file is parsed and the whole representation of file storage is built in memory as a hierarchy of file nodes (see FileNode)
// 2. Read the data you are interested in. Use FileStorage::operator [](), FileNode::operator [](), and/or FileNodeIterator.
//   a) first method: use (type) operator on FileNode;
//   b) second method: use FileNode::operator >>;
// 3. Close the storage using FileStorage::release().

#include "opencv2/opencv.hpp"
#include <iostream>
#include <vector>
using namespace cv;
using namespace std;
#define Read_FileStorage 1
#define Write_FileStorage 0

/*1. static function as a plain c function:
static functions are functions that are only visible to other functions in the same file.

2. static function as a class method:
For class methods, static obviously means that this method can be called on the class itself, no instance of that class necessary.
*/

void fileStorageHelp();

void writeDataToFileStorage(vector<vector<vector<int>>> & v,
	string & filename, // e.g. filename = "test.yml"
	string & fileNodeName // e.g. "imageIdx"
	);



void readDataFromFileStorage(vector<vector<vector<int>>> & v,
	string & filename, // e.g. filename = "test.yml"
	string & fileNodeName // e.g. "imageIdx"
	);


void read_Write_FileStorage(vector<vector<vector<int>>> & v,
	string & filename, // e.g. filename = "test.yml"
	string & filenodename,
	int flag  // #define Read_FileStorage 1; #define Write_FileStorage 0
	);

#endif