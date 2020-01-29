#include "convolution_std.h"
#include <ctime> // clock_t t = clock();
namespace STD_Convolution{

	void init_workspace(Std_Workspace & ws, Std_Convolution_Mode mode, int h_src, int w_src, int h_kernel, int w_kernel, bool verbose)
	{
		ws.h_src = h_src;
		ws.w_src = w_src;
		ws.h_kernel = h_kernel;
		ws.w_kernel = w_kernel;
		ws.mode = mode;
		ws.showSeconds = verbose;

		switch (ws.mode){
		case 0:
			ws.s_mode = "LINEAR_CONVOLUTION_FULL"; 
			break;
		
		case 1:
			ws.s_mode = "LINEAR_CONVOLUTION_SAME";
			break;
		
		case 2:
			ws.s_mode = "LINEAR_CONVOLUTION_VALID";
			break;
		
		case 3:
			ws.s_mode = "LINEAR_CONVOLUTION_VALIDSAME";
			break;
		
		case 4:
			ws.s_mode = "LINEAR_CORRELATION_FULL";
			break;
		
		case 5:
			ws.s_mode = "LINEAR_CORRELATION_SAME";
			break;
		
		case 6:
			ws.s_mode = "LINEAR_CORRELATION_VALID";
			break;
		
		case 7:
			ws.s_mode = "LINEAR_CORRELATION_VALIDSAME";
			break;
		
		case 8:
			ws.s_mode = "CIRCULAR_CONVOLUTION_FULL";
			break;
		
		case 9:
			ws.s_mode = "CIRCULAR_CONVOLUTION_SAME";
			break;
		
		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
            printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
		}


		switch (mode){
		case LINEAR_CONVOLUTION_FULL:
		case LINEAR_CORRELATION_FULL:
			ws.h_dst = ws.h_src + ws.h_kernel - 1;
			ws.w_dst = ws.w_src + ws.w_kernel - 1;
			ws.offset_h_dst = 0;
			ws.offset_w_dst = 0;
			break;
		case LINEAR_CONVOLUTION_SAME:
		case LINEAR_CORRELATION_SAME:
			ws.h_dst = ws.h_src;
			ws.w_dst = ws.w_src;
			ws.offset_h_dst = int(ws.h_kernel /2.0);
			ws.offset_w_dst = int(ws.w_kernel /2.0);
			break;

		case LINEAR_CONVOLUTION_VALID:
		case LINEAR_CORRELATION_VALID:
			if (ws.h_kernel > ws.h_src || ws.w_kernel > ws.w_src){
				//printf("Warning : The 'valid' convolution results in an empty matrix\n");
				ws.h_dst = 0;
				ws.w_dst = 0;
				ws.offset_h_dst = ws.h_kernel - 1;
				ws.offset_w_dst = ws.w_kernel - 1;
			}
			else{
				ws.h_dst = ws.h_src - ws.h_kernel + 1;
				ws.w_dst = ws.w_src - ws.w_kernel + 1;
				ws.offset_h_dst = ws.h_kernel - 1;
				ws.offset_w_dst = ws.w_kernel - 1;
			}
			break;

		case LINEAR_CONVOLUTION_VALIDSAME: // NON PADDING
		case LINEAR_CORRELATION_VALIDSAME: // NON PADDING
			if (ws.h_kernel > ws.h_src || ws.w_kernel > ws.w_src){
				//printf("Warning : The 'valid' convolution results in an empty matrix\n");
				ws.h_dst = 0;
				ws.w_dst = 0;
				ws.offset_h_dst = ws.h_kernel - 1;
				ws.offset_w_dst = ws.w_kernel - 1;
			}
			else{
				ws.h_dst = h_src;
				ws.w_dst = w_src;
				ws.offset_h_dst = ws.h_kernel - 1;
				ws.offset_w_dst = ws.w_kernel - 1;
			}
			break;

		case CIRCULAR_CONVOLUTION_SAME:
			ws.h_dst = ws.h_src;
			ws.w_dst = ws.w_src;
			
			break;
		case CIRCULAR_CONVOLUTION_FULL:
			ws.h_dst = ws.h_src + ws.h_kernel - 1;
			ws.w_dst = ws.w_src + ws.w_kernel - 1;
			break;

		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
		}

		if (ws.h_dst > 0 && ws.w_dst > 0)
			ws.dst = new double[ws.h_dst * ws.w_dst];
		//else
		//  printf("Warning : The result is an empty matrix !\n");

	}

	void clear_workspace(Std_Workspace & ws) // to delete output array
	{
		if (ws.h_dst > 0 && ws.w_dst > 0)
			delete[] ws.dst;
	}

	void update_workspace(Std_Workspace & ws, Std_Convolution_Mode mode, int h_src, int w_src, int h_kernel, int w_kernel, bool verbose)
	{
		clear_workspace(ws);
		init_workspace(ws, mode, h_src, w_src, h_kernel, w_kernel, verbose );
	}



	void convolve(Std_Workspace &ws, double * src, double * kernel)
	{
		clock_t t = clock();
		double temp;
		int row_Idx_dst, col_Idx_dst; // for output image
		int row_Idx_src, col_Idx_src; // for input image 
		// lower bound and upper bound for input image
		int low_RowIdxSrc, high_RowIdxSrc; // for row indices
		int low_ColIdxSrc, high_ColIdxSrc; // for column indices
		int circular_row_Idx_src, circular_col_Idx_src; // circular indices of input images for circular convolution

		if (ws.h_dst <= 0 || ws.w_dst <= 0)
			return;

#define TRY_SWITCH_CASE_WITH_FUNCTION
#ifdef TRY_SWITCH_CASE_WITH_FUNCTION
		switch (ws.mode){
			// for correlation, we have to flip the kernel
		case LINEAR_CORRELATION_FULL: // 4		
		case LINEAR_CORRELATION_SAME: // 5		
		case LINEAR_CORRELATION_VALID: // 6		
		case LINEAR_CORRELATION_VALIDSAME: // 7
			// first flip the kernel
			std::reverse(kernel, kernel + ws.h_kernel*ws.w_kernel);
			// then do the same thing to the src and the flipped kernel as the case of convolution does.
			general_convolve(ws,src, kernel);
			break;

		default:// for convolution, we DO NOT flip the kernel
			// so just directly do the convolution.
			general_convolve(ws, src, kernel);
		}
#endif

#define TRY_SWITCH_CASE_NO_FUNCTION
#ifndef TRY_SWITCH_CASE_NO_FUNCTION

		switch (ws.mode)
		{
		case LINEAR_CONVOLUTION_FULL: // Full linear convolution of size N + M -1

			// see http://stackoverflow.com/questions/484462/difference-between-i-and-i-in-a-loop
			// In C++ if you're using STL, then you may be using for loops with iterators. 
			// These mainly have overridden ++ operators, so sticking to pre-increment is a good idea. 

			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst) { // loop for output image (i.e., convolution result) height
				// /* by CCJ
				//  * see http://stackoverflow.com/questions/7035023/stdmax-expected-an-identifier
				//  * Answer 2 (CCJ recommends it)
				//  * e.g.
				//  * int x = 2, y = 4;
				//  * std::max(x, y) // to cause a problem, like "std::max expected an identifier", due to the conflict between std::max and windows.h
				//  * since windows.h has  " #define max(a,b)   (((a) > (b)) ? (a) : (b))" 
				//  * std::max<int>(x,y) // is OK
				//  * (std::max)(x, y) // is OK
				//  */
				low_RowIdxSrc = std::max<int>(0, row_Idx_dst - ws.h_kernel + 1); // need to explicitly invoke the template via std::max<int> 
				high_RowIdxSrc = std::min<int>(ws.h_src - 1, row_Idx_dst);
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst)   // loop for output image (i.e., convolution result) width
				{
					low_ColIdxSrc = (std::max)(0, col_Idx_dst - ws.w_kernel + 1); // (std::max)(...) is OK.
					high_ColIdxSrc = (std::min)(ws.w_src - 1, col_Idx_dst);

					temp = 0.0;
					for (row_Idx_src = low_RowIdxSrc; row_Idx_src <= high_RowIdxSrc; ++row_Idx_src){  // loop for input image's height
						for (col_Idx_src = low_ColIdxSrc; col_Idx_src <= high_ColIdxSrc; ++col_Idx_src){  // loop for input image's width
							temp += src[row_Idx_src*ws.w_src + col_Idx_src] * kernel[(row_Idx_dst - row_Idx_src)*ws.w_kernel + (col_Idx_dst - col_Idx_src)];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		case LINEAR_CORRELATION_FULL: // Full linear correlation of size N + M -1
			// first flip the kernel
			std::reverse(kernel, kernel + ws.h_kernel*ws.w_kernel);
			// then do the same thing to the src and the flipped kernel as the case of LINEAR_CONVOLUTION_FULL does.
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst) { // loop for output image (i.e., convolution result) height
				low_RowIdxSrc = std::max<int>(0, row_Idx_dst - ws.h_kernel + 1); // need to explicitly invoke the template via std::max<int> 
				high_RowIdxSrc = std::min<int>(ws.h_src - 1, row_Idx_dst);
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst)   // loop for output image (i.e., convolution result) width
				{
					low_ColIdxSrc = (std::max)(0, col_Idx_dst - ws.w_kernel + 1); // (std::max)(...) is OK.
					high_ColIdxSrc = (std::min)(ws.w_src - 1, col_Idx_dst);

					temp = 0.0;
					for (row_Idx_src = low_RowIdxSrc; row_Idx_src <= high_RowIdxSrc; ++row_Idx_src){  // loop for input image's height
						for (col_Idx_src = low_ColIdxSrc; col_Idx_src <= high_ColIdxSrc; ++col_Idx_src){  // loop for input image's width
							temp += src[row_Idx_src*ws.w_src + col_Idx_src] * kernel[(row_Idx_dst - row_Idx_src)*ws.w_kernel + (col_Idx_dst - col_Idx_src)];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		

		case LINEAR_CONVOLUTION_SAME:
			// Same linear convolution, of size N
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				low_RowIdxSrc = std::max<int>(0, row_Idx_dst - int((ws.h_kernel - 1.0) / 2.0));
				high_RowIdxSrc = std::min<int>(ws.h_src - 1, row_Idx_dst + int(ws.h_kernel / 2.0));
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst) { // loop for output image's (i.e., convolution result) width
					low_ColIdxSrc = std::max<int>(0, col_Idx_dst - int((ws.w_kernel - 1.0) / 2.0));
					high_ColIdxSrc = std::min<int>(ws.w_src - 1, col_Idx_dst + int(ws.w_kernel / 2.0));
					temp = 0.0;
					for (row_Idx_src = low_RowIdxSrc; row_Idx_src <= high_RowIdxSrc; ++row_Idx_src){ // loop for input image's height
						for (col_Idx_src = low_ColIdxSrc; col_Idx_src <= high_ColIdxSrc; ++col_Idx_src){ // loop for input image's width
							temp += src[row_Idx_src*ws.w_src + col_Idx_src] * kernel[(row_Idx_dst - row_Idx_src + int(ws.h_kernel / 2.0))*ws.w_kernel + (col_Idx_dst - col_Idx_src + int(ws.w_kernel / 2.0))];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		case LINEAR_CORRELATION_SAME: // SAME linear correlation of size N
			// first flip the kernel
			std::reverse(kernel, kernel + ws.h_kernel*ws.w_kernel);
			// then do the same thing to the src and the flipped kernel as the case of LINEAR_CORRELATION_SAME does.
			// Same as the Same-linear convolution, of size N
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				low_RowIdxSrc = std::max<int>(0, row_Idx_dst - int((ws.h_kernel - 1.0) / 2.0));
				high_RowIdxSrc = std::min<int>(ws.h_src - 1, row_Idx_dst + int(ws.h_kernel / 2.0));
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst) { // loop for output image's (i.e., convolution result) width
					low_ColIdxSrc = std::max<int>(0, col_Idx_dst - int((ws.w_kernel - 1.0) / 2.0));
					high_ColIdxSrc = std::min<int>(ws.w_src - 1, col_Idx_dst + int(ws.w_kernel / 2.0));
					temp = 0.0;
					for (row_Idx_src = low_RowIdxSrc; row_Idx_src <= high_RowIdxSrc; ++row_Idx_src){ // loop for input image's height
						for (col_Idx_src = low_ColIdxSrc; col_Idx_src <= high_ColIdxSrc; ++col_Idx_src){ // loop for input image's width
							temp += src[row_Idx_src*ws.w_src + col_Idx_src] * kernel[(row_Idx_dst - row_Idx_src + int(ws.h_kernel / 2.0))*ws.w_kernel + (col_Idx_dst - col_Idx_src + int(ws.w_kernel / 2.0))];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		case LINEAR_CONVOLUTION_VALID:
			// Valid linear convolution, of size (N - M + 1)
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst){ // loop for output image's (i.e., convolution result) width
					temp = 0.0;
					for (row_Idx_src = row_Idx_dst; row_Idx_src <= row_Idx_dst + ws.h_kernel - 1; ++row_Idx_src){ // loop for input image's height
						for (col_Idx_src = col_Idx_dst; col_Idx_src <= col_Idx_dst + ws.w_kernel - 1; ++col_Idx_src){ // loop for input image's width
							temp += src[row_Idx_src*ws.w_src + col_Idx_src] * kernel[(row_Idx_dst + ws.h_kernel - 1 - row_Idx_src)*ws.w_kernel + (col_Idx_dst + ws.w_kernel - 1 - col_Idx_src)];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;
		case LINEAR_CORRELATION_VALID:
			// first flip the kernel
			std::reverse(kernel, kernel + ws.h_kernel*ws.w_kernel);
			// then do the same thing to the src and the flipped kernel as the case of LINEAR_CONVOLUTION_VALID does.
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst){ // loop for output image's (i.e., convolution result) width
					temp = 0.0;
					for (row_Idx_src = row_Idx_dst; row_Idx_src <= row_Idx_dst + ws.h_kernel - 1; ++row_Idx_src){ // loop for input image's height
						for (col_Idx_src = col_Idx_dst; col_Idx_src <= col_Idx_dst + ws.w_kernel - 1; ++col_Idx_src){ // loop for input image's width
							temp += src[row_Idx_src*ws.w_src + col_Idx_src] * kernel[(row_Idx_dst + ws.h_kernel - 1 - row_Idx_src)*ws.w_kernel + (col_Idx_dst + ws.w_kernel - 1 - col_Idx_src)];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		

		case LINEAR_CONVOLUTION_VALIDSAME: // NOTE: ONLY for the VALIDSAME case, we do padding to input signal via circular shifting.
			// Valid&Same linear convolution, of size (N )
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst){ // loop for output image's (i.e., convolution result) width
					temp = 0.0;
					for (row_Idx_src = row_Idx_dst; row_Idx_src <= row_Idx_dst + ws.h_kernel - 1; ++row_Idx_src){ // loop for input image's height
						for (col_Idx_src = col_Idx_dst; col_Idx_src <= col_Idx_dst + ws.w_kernel - 1; ++col_Idx_src){ // loop for input image's width
							// here we do circular-shifting to src, with modulo N, when its indices is larger than N.
							temp += src[(row_Idx_src % ws.h_src)*ws.w_src + (col_Idx_src % ws.w_src)] * kernel[(row_Idx_dst + ws.h_kernel - 1 - row_Idx_src)*ws.w_kernel + (col_Idx_dst + ws.w_kernel - 1 - col_Idx_src)];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;
		
		case LINEAR_CORRELATION_VALIDSAME: // NOTE: ONLY for the VALIDSAME case, we do padding to input signal via circular shifting.
			// first flip the kernel
			std::reverse(kernel, kernel + ws.h_kernel*ws.w_kernel);
			// then do the same thing to the src and the flipped kernel as the case of LINEAR_CONVOLUTION_VALIDSAME does.
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst){ // loop for output image's (i.e., convolution result) width
					temp = 0.0;
					for (row_Idx_src = row_Idx_dst; row_Idx_src <= row_Idx_dst + ws.h_kernel - 1; ++row_Idx_src){ // loop for input image's height
						for (col_Idx_src = col_Idx_dst; col_Idx_src <= col_Idx_dst + ws.w_kernel - 1; ++col_Idx_src){ // loop for input image's width
							// here we do circular-shifting to src, with modulo N, when its indices is larger than N.
							temp += src[(row_Idx_src % ws.h_src)*ws.w_src + (col_Idx_src % ws.w_src)] * kernel[(row_Idx_dst + ws.h_kernel - 1 - row_Idx_src)*ws.w_kernel + (col_Idx_dst + ws.w_kernel - 1 - col_Idx_src)];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		case CIRCULAR_CONVOLUTION_SAME:
		case CIRCULAR_CONVOLUTION_FULL:
			// We suppose the filter has a size at most the size of the image
			// i.e., size of kernel <= size of src
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst){ // loop for output image's (i.e., convolution result) width
					temp = 0.0;

					// We browse the kernel
					// here I use (*) to present circular convolution
					// Since, (input (*) kernel ) [k] = (kernel (*) input)[k], 
					// (input (*) kernel ) [k], for any k belonging to [0, N-1], with N being the length of input signal or image,
					// equals (kernel (*) input)[k], for any k belonging to [0, M-1], with M being the length of kernel;
					for (row_Idx_src = 0; row_Idx_src < ws.h_kernel; ++row_Idx_src){ // kernel's height
						circular_row_Idx_src = (row_Idx_dst - row_Idx_src) % ws.h_dst;  // this is the formula after we change the order of input and kernel signals 
						if (circular_row_Idx_src < 0)
							circular_row_Idx_src += ws.h_dst;

						for (col_Idx_src = 0; col_Idx_src < ws.w_kernel; ++col_Idx_src){ // kernel's width
							circular_col_Idx_src = (col_Idx_dst - col_Idx_src) % ws.w_dst; // this is the formula after we change the order of input and kernel signals
							if (circular_col_Idx_src < 0)
								circular_col_Idx_src += ws.w_dst;


							// /* for circular convolution, there also exist two cases, i.e., CIRCULAR_FULL and CIRCULAR_SAME
							//  * in the following if-sentence we have two condition "circular_row_Idx_src < ws.h_src" and "circular_col_Idx_src < ws.w_src" in case of CIRCULAR_FULL
							//  * they are used to control the indices for the CIRCULAR_FULL case;
							//  * for CIRCULAR_SAME, the indices will not exceed the bounds of input image "src"
							//  */
							if (circular_row_Idx_src >= 0 && circular_row_Idx_src < ws.h_src  && circular_col_Idx_src >= 0 && circular_col_Idx_src < ws.w_src)
								temp += kernel[row_Idx_src*ws.w_kernel + col_Idx_src] * src[circular_row_Idx_src*ws.w_src + circular_col_Idx_src] ;
						}
					}
					ws.dst[row_Idx_dst*ws.w_dst + col_Idx_dst] = temp;
				}
			}

			break;
		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			break;
		} // end of switch-case
#endif
		if (ws.showSeconds){
			t = clock() - t;
			printf("STD_Convolution::convolve function took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
		}
}

	
// find a general function to calculate the different kinds of convolution or correlation
void general_convolve(Std_Workspace & ws, double * src, double * kernel){
	    clock_t t = clock();
		double temp;
		int row_Idx_dst, col_Idx_dst; // for output image
		int row_Idx_src, col_Idx_src; // for input image 
		
		// lower bound and upper bound for input image
		int low_RowIdxSrc, high_RowIdxSrc; // for row indices
		int low_ColIdxSrc, high_ColIdxSrc; // for column indices
		
		int circular_row_Idx_src, circular_col_Idx_src; // circular indices of input images for circular convolution

		switch (ws.mode){
		case LINEAR_CORRELATION_VALIDSAME:
		case LINEAR_CONVOLUTION_VALIDSAME:
			// NOTE: ONLY for the VALIDSAME case, we do padding to input signal via circular shifting.
			// Valid&Same linear convolution, of size (N )
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst){ // loop for output image's (i.e., convolution result) width
					temp = 0.0;
					for (row_Idx_src = row_Idx_dst; row_Idx_src <= row_Idx_dst + ws.h_kernel - 1; ++row_Idx_src){ // loop for input image's height
						for (col_Idx_src = col_Idx_dst; col_Idx_src <= col_Idx_dst + ws.w_kernel - 1; ++col_Idx_src){ // loop for input image's width
							// here we do circular-shifting to src, with modulo N, when its indices is larger than N.
							temp += src[(row_Idx_src % ws.h_src)*ws.w_src + (col_Idx_src % ws.w_src)] * kernel[(row_Idx_dst + ws.h_kernel - 1 - row_Idx_src)*ws.w_kernel + (col_Idx_dst + ws.w_kernel - 1 - col_Idx_src)];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		case CIRCULAR_CONVOLUTION_SAME:
		case CIRCULAR_CONVOLUTION_FULL:
			// We suppose the filter has a size at most the size of the image
			// i.e., size of kernel <= size of src
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst){ // loop for output image's (i.e., convolution result) height
				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst){ // loop for output image's (i.e., convolution result) width
					temp = 0.0;

					// We browse the kernel
					// here I use (*) to present circular convolution
					// Since, (input (*) kernel ) [k] = (kernel (*) input)[k], 
					// (input (*) kernel ) [k], for any k belonging to [0, N-1], with N being the length of input signal or image,
					// equals (kernel (*) input)[k], for any k belonging to [0, M-1], with M being the length of kernel;
					for (row_Idx_src = 0; row_Idx_src < ws.h_kernel; ++row_Idx_src){ // kernel's height
						circular_row_Idx_src = (row_Idx_dst - row_Idx_src) % ws.h_dst;  // this is the formula after we change the order of input and kernel signals 
						if (circular_row_Idx_src < 0)
							circular_row_Idx_src += ws.h_dst;

						for (col_Idx_src = 0; col_Idx_src < ws.w_kernel; ++col_Idx_src){ // kernel's width
							circular_col_Idx_src = (col_Idx_dst - col_Idx_src) % ws.w_dst; // this is the formula after we change the order of input and kernel signals
							if (circular_col_Idx_src < 0)
								circular_col_Idx_src += ws.w_dst;


							// /* for circular convolution, there also exist two cases, i.e., CIRCULAR_FULL and CIRCULAR_SAME
							//  * in the following if-sentence we have two condition "circular_row_Idx_src < ws.h_src" and "circular_col_Idx_src < ws.w_src" in case of CIRCULAR_FULL
							//  * they are used to control the indices for the CIRCULAR_FULL case;
							//  * for CIRCULAR_SAME, the indices will not exceed the bounds of input image "src"
							//  */
							if (circular_row_Idx_src >= 0 && circular_row_Idx_src < ws.h_src  && circular_col_Idx_src >= 0 && circular_col_Idx_src < ws.w_src)
								temp += kernel[row_Idx_src*ws.w_kernel + col_Idx_src] * src[circular_row_Idx_src*ws.w_src + circular_col_Idx_src];
						}
					}
					ws.dst[row_Idx_dst*ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		case LINEAR_CONVOLUTION_FULL: // 0
		case LINEAR_CONVOLUTION_SAME: // 1
		case LINEAR_CONVOLUTION_VALID: // 2
		case LINEAR_CORRELATION_FULL: // 4		
		case LINEAR_CORRELATION_SAME: // 5		
		case LINEAR_CORRELATION_VALID: // 6
			// otherwise, then, without padding to input signal via circular shifting.
			for (row_Idx_dst = 0; row_Idx_dst < ws.h_dst; ++row_Idx_dst) { // loop for output image (i.e., convolution result) height
				// /* by CCJ
				//  * see http://stackoverflow.com/questions/7035023/stdmax-expected-an-identifier
				//  * Answer 2 (CCJ recommends it)
				//  * e.g.
				//  * int x = 2, y = 4;
				//  * std::max(x, y) // to cause a problem, like "std::max expected an identifier", due to the conflict between std::max and windows.h
				//  * since windows.h has  " #define max(a,b)   (((a) > (b)) ? (a) : (b))" 
				//  * std::max<int>(x,y) // is OK
				//  * (std::max)(x, y) // is OK
				//  */

				low_RowIdxSrc = std::max<int>(0, row_Idx_dst + ws.offset_h_dst - ws.h_kernel + 1); // need to explicitly invoke the template via std::max<int> 
				high_RowIdxSrc = std::min<int>(ws.h_src - 1, row_Idx_dst + ws.offset_h_dst);

				for (col_Idx_dst = 0; col_Idx_dst < ws.w_dst; ++col_Idx_dst)   // loop for output image (i.e., convolution result) width
				{
					low_ColIdxSrc = (std::max)(0, col_Idx_dst + ws.offset_w_dst - ws.w_kernel + 1); // (std::max)(...) is OK.
					high_ColIdxSrc = (std::min)(ws.w_src - 1, col_Idx_dst + ws.offset_w_dst);

					temp = 0.0;
					for (row_Idx_src = low_RowIdxSrc; row_Idx_src <= high_RowIdxSrc; ++row_Idx_src){  // loop for input image's height
						for (col_Idx_src = low_ColIdxSrc; col_Idx_src <= high_ColIdxSrc; ++col_Idx_src){  // loop for input image's width
							temp += src[row_Idx_src*ws.w_src + col_Idx_src] * kernel[(row_Idx_dst + ws.offset_h_dst - row_Idx_src)*ws.w_kernel + (col_Idx_dst + ws.offset_w_dst - col_Idx_src)];
						}
					}
					ws.dst[row_Idx_dst * ws.w_dst + col_Idx_dst] = temp;
				}
			}
			break;

		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			break;
		}
		if (ws.showSeconds){
			t = clock() - t;
		    printf("STD_Convolution::general_convolve function took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
		}
	}

	
	void saveConvolve(Std_Workspace & ws, std::string & filename){
		std::ofstream fout(filename, std::ios::app);
		// /*
		//  * ios::in  --- Open for input operations;
		//  * ios::out --- Open for output operations;
		//  * ios::binary ---- Open in binary mode.
		//  * ios::ate ---- Set the initial position at the end of the file. If this flag is not set, the initial position is the beginning of the file.
		//  * ios::app --- All output operations are performed at the end of the file, appending the content to the current content of the file.
		//  * ios::trunc ---- If the file is opened for output operations and it already existed, its previous content is deleted and replaced by the new one.
		//  */

		if (!fout){
			std::cout << "File Not Opened" << std::endl;
		}
		fout << "Input Signal : " << ws.h_src << " by " << ws.w_src << std::endl
			<< "Standard Convolution/Correlation mode is : " << ws.s_mode << std::endl
			<< "Standard Convolution/Correlation 2-D result : " << std::endl;
		for (int i = 0; i < ws.h_dst; ++i){
			for (int j = 0; j < ws.w_dst; ++j){
				fout << ws.dst[i*ws.w_dst + j] << ", ";
			}
			fout << std::endl;
		}
		fout << std::endl;
		fout.close();
	}
}