/*---------------------------------------------------------------------------*/
// Image Compression Toolbox v1.2
// written by
// Satish Kumar S
// satishkumr@lycos.com
//
// Copyright 1999 Satish Kumar S
//
// Permission is granted to use this software for research purposes as
// long as this notice stays attached to this software.
/*---------------------------------------------------------------------------*/

#include "bitstream.h"

//-----------------------------------------------------------
// the Bit Stream part.
//-----------------------------------------------------------
BitStream::BitStream(FILE *givenfp){
	buf = 0;
	nBits = 0;
	numBitsWritten = 0;
	numBitsRead = 0;
	fp = givenfp;
}

BitStream::~BitStream()
{
	buf = 0;
	nBits = 0;
	numBitsWritten = 0; 
	numBitsRead = 0;
	fp = NULL;
}

int BitStream::SetFilePointer(FILE * newfp){
	if (newfp)
		fp = newfp;
	else
		return 1;
	return 0;
}

void BitStream::Flush(){

	// negative or positive. Anything that's not a 0 is a true value in if statement
	if (nBits){ // nonzero value means TRUE
		buf = (buf << (8 - nBits));
		fwrite(&buf, 1, 1, fp);
		nBits = 0;
		buf = 0;
	}
}

void BitStream::clear(){
	buf = 0;
	nBits = 0;
	numBitsWritten = 0;
	numBitsRead = 0;
	fp = NULL;
}

// write only one bitValue ( 0 or 1) to buffer buf.
// this function will be called by BitStream::WriteBits(...) to finish more complicated task.
void BitStream::WriteBit(unsigned int bitValue // bitValue == 0 or 1
	){

	// first left-shift buf by 1 bit, then add the bitValue (i.e., 0 or 1)
	buf = (buf << 1) + bitValue;
	nBits++;

	if (nBits == 8){ // if 8 bits have been achieved, then write them into the file pointed by fp
		fwrite(&buf, 1, 1, fp);
		nBits = 0;
		buf = 0;
	}
	numBitsWritten++;
}

// this function will call BitStream::WriteBit(...) to finish writing the lowest n bits of the integer of "bits"
// into the buffer "buf"
void BitStream::WriteBits(unsigned int bits, int n){
	unsigned int mask = (1 << (n - 1));
	while (n--){
		// The bitwise AND operator : &
		// relational operator : > (less than)
		// 关系运算符的值只能是0或1。
		// 关系运算符的值为真时，结果值都为1。
		// 关系运算符的值为假时，结果值都为0。
		// that is, the value of (bits & mask) > 0, will be 1 (true), or 0 (false).
		WriteBit((bits & mask) > 0);
		bits = bits << 1;
	}
}

// read the highest 1 bit of the value in buf
// 即，一次只读一个位，并且是最高位。
// this function will be called by BitStream::ReadBits(int n), to finish n-bit reading from file.
unsigned int BitStream::ReadBit(){
	// if nBits == 0 , that means we have read all the 8-bit of the unsigned char value in buf
	// then, we have to read a new number to buf from file, and set nBits == 8 for next 8-bit reading cycle. 
	if (nBits == 0){
		nBits = 8;
		// function : size_t fread ( void * ptr, size_t size, size_t count, FILE * stream );
		// Reads an array of count elements, each one with a size of size bytes, from the stream and stores them in the block of memory specified by ptr.
		// The position indicator of the stream is advanced by the total amount of bytes read.
		// The total amount of bytes read if successful is (size*count).
		fread(&buf, 1, 1, fp);
	}

	// 128 == 0x80 == 1000 0000 (binary number)
	// to get the highest bit in buf, and then shift it from the highest bit to the lowest bit.
	int bit = (buf & 128) >> 7;

	// left shift buf by 1 bit
	buf = buf << 1;
	nBits--;
	numBitsRead++;
	return bit;
}


// this function will call BitStream::ReadBit(), 
// to extract the highest n-bit of the value in buf (but the value in buf has been read from file)
// returns an unsigned int number
unsigned int BitStream::ReadBits(int n){
	int bits = 0;
	while (n--)
		bits = (bits << 1) + ReadBit();
	return bits;
}