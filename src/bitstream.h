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
#ifndef _BIT_STREAM_H
#define _BIT_STREAM_H

#include <stdio.h>

//--------- Bit Stream part.
class BitStream{

private:
	unsigned char buf;
	int nBits;
	FILE *fp;
public:
	int  numBitsWritten;
	int  numBitsRead;
	BitStream(FILE *givenfp);
	~BitStream();
	void WriteBit(unsigned int bit);
	void WriteBits(unsigned int bits, int n);
	void Flush();
	unsigned int ReadBit();
	unsigned int ReadBits(int n);
	int SetFilePointer(FILE * newfp);
	void BitStream::clear();
};

#endif