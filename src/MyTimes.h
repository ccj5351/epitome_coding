// From http://stackoverflow.com/questions/10905892/equivalent-of-gettimeday-for-windows

#ifndef _MY_TIMES_H_
#define _MY_TIMES_H_

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <stdint.h> // portable: uint64_t   MSVC: __int64 


// MSVC defines this in winsock2.h!?
struct timeval{
	long tv_sec;
	long tv_usec;
};

int gettimeofday( timeval * tp, struct timezone * tzp);


#endif