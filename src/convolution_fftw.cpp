#include "convolution_fftw.h"

namespace FFTW_Convolution{

	void init_workspace(FFT_Workspace & ws, FFT_Convolution_Mode mode, int h_src, int w_src, int h_kernel, int w_kernel, bool verbose)
	{
		ws.showSeconds = verbose;
		ws.h_src = h_src;
		ws.w_src = w_src;
		ws.h_kernel = h_kernel;
		ws.w_kernel = w_kernel;
		ws.mode = mode;
		
		// If your program performs many transforms of the same size and initialization time is not important, use FFTW_MEASURE; 
		// otherwise use FFTW_ESTIMATE.
		ws.fftw_flag = FFTW_MEASURE;
		// ws.fftw_flag = FFTW_MEASURE;

		// FFTW_FACTORS = { 13, 11, 7, 5, 3, 2, 0 }; // end with zero to detect the end of the array
		ws.FFTW_FACTORS = new int[7];
		ws.FFTW_FACTORS[0] = 13;
		ws.FFTW_FACTORS[1] = 11;
		ws.FFTW_FACTORS[2] = 7;
		ws.FFTW_FACTORS[3] = 5;
		ws.FFTW_FACTORS[4] = 3;
		ws.FFTW_FACTORS[5] = 2;
		ws.FFTW_FACTORS[6] = 0;

		switch (ws.mode){
		case 0:
			ws.s_mode = "LINEAR_CONVOLUTION_FULL"; break;
		case 10:
			ws.s_mode = "LINEAR_CONVOLUTION_FULL_UNPADDED"; break;
		case 1:
			ws.s_mode = "LINEAR_CONVOLUTION_SAME";
			break;
		case 11:
			ws.s_mode = "LINEAR_CONVOLUTION_SAME_UNPADDED";
			break;
		case 2:
			ws.s_mode = "LINEAR_CONVOLUTION_VALID";
			break;
		case 12:
			ws.s_mode = "LINEAR_CONVOLUTION_VALID_UNPADDED";
			break;
		case 3:
			ws.s_mode = "LINEAR_CONVOLUTION_VALIDSAME";
			break;
		case 13:
			ws.s_mode = "LINEAR_CONVOLUTION_VALIDSAME_UNPADDED";
			break;
		case 4:
			ws.s_mode = "LINEAR_CORRELATION_FULL";
			break;
		case 14:
			ws.s_mode = "LINEAR_CORRELATION_FULL_UNPADDED";
			break;
		case 5:
			ws.s_mode = "LINEAR_CORRELATION_SAME";
			break;
		case 15:
			ws.s_mode = "LINEAR_CORRELATION_SAME_UNPADDED";
			break;
		case 6:
			ws.s_mode = "LINEAR_CORRELATION_VALID";
			break;
		case 16:
			ws.s_mode = "LINEAR_CORRELATION_VALID_UNPADDED";
			break;
		case 7:
			ws.s_mode = "LINEAR_CORRELATION_VALIDSAME";
			break;
		case 17:
			ws.s_mode = "LINEAR_CORRELATION_VALIDSAME_UNPADDED";
			break;
		case 8:
			ws.s_mode = "CIRCULAR_CONVOLUTION_FULL";
			break; 
		case 18:
				ws.s_mode = "CIRCULAR_CONVOLUTION_FULL_UNPADDED";
				break;
		case 9:
			ws.s_mode = "CIRCULAR_CONVOLUTION_SAME";
			break;
		case 19:
			ws.s_mode = "CIRCULAR_CONVOLUTION_SAME_UNPADDED";
			break;
		case 20:
			ws.s_mode = "LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP";
			break;
		case 21:
			ws.s_mode = "LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP";
			break;
		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_SAME_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALID_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME_UNPADDED \n");

			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_FULL_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_SAME_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALID_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_UNPADDED \n");

			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME_UNPADDED \n");

			printf("   - LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP \n");		
		}
		
		
		switch (ws.mode){
			// Here if no explicit "UNPADDED" occurs, we always mean PADDING as default.
			// i.e., here LINEAR_FULL == LINEAR_FULL_PADDED
		
			// *********************************************
			// Full Linear convolution or correlation
			// *********************************************
		case LINEAR_CONVOLUTION_FULL: // PADDING
		case LINEAR_CORRELATION_FULL: // PADDING
			ws.h_fftw = find_closest_factor(h_src + h_kernel - 1, ws.FFTW_FACTORS);
			ws.w_fftw = find_closest_factor(w_src + w_kernel - 1, ws.FFTW_FACTORS);
			ws.h_dst = h_src + h_kernel - 1;
			ws.w_dst = w_src + w_kernel - 1;
			ws.IsPeriodicPadding2Src = false;
		//	ws.IsCircularlyWrapForValidSame = false;
			break;
			
		case LINEAR_CONVOLUTION_FULL_UNPADDED: // NON PADDING
		case LINEAR_CORRELATION_FULL_UNPADDED: // NON PADDING
			ws.h_dst = h_src + h_kernel - 1;
			ws.w_dst = w_src + w_kernel - 1;
			ws.h_fftw = ws.h_dst;
			ws.w_fftw = ws.w_dst;
			ws.IsPeriodicPadding2Src = false;
		//	ws.IsCircularlyWrapForValidSame = false;
			break;


			// *********************************************
			// Same Linear convolution or correlation
			// *********************************************
		case LINEAR_CONVOLUTION_SAME_UNPADDED: // NO PADDING
		case LINEAR_CORRELATION_SAME_UNPADDED: // NO PADDING
			ws.h_dst = h_src;
			ws.w_dst = w_src;
			ws.h_fftw = h_src + int(h_kernel / 2.0);
			ws.w_fftw = w_src + int(w_kernel / 2.0);
			ws.IsPeriodicPadding2Src = false;
		//	ws.IsCircularlyWrapForValidSame = false;
			break;

		case LINEAR_CONVOLUTION_SAME: // PADDING
		case LINEAR_CORRELATION_SAME: // PADDING
			ws.h_fftw = find_closest_factor(h_src + int(h_kernel / 2.0), ws.FFTW_FACTORS);
			ws.w_fftw = find_closest_factor(w_src + int(w_kernel / 2.0), ws.FFTW_FACTORS);
			ws.h_dst = h_src;
			ws.w_dst = w_src;
			ws.IsPeriodicPadding2Src = false;
		//	ws.IsCircularlyWrapForValidSame = false;
			break;

			// *********************************************
			// Valid Linear convolution or correlation
			// *********************************************
		case LINEAR_CONVOLUTION_VALID: // PADDING
		case LINEAR_CORRELATION_VALID: // PADDING
			if (ws.h_kernel > ws.h_src || ws.w_kernel > ws.w_src){
				//printf("Warning : The 'valid' convolution results in an empty matrix\n");
				ws.h_fftw = 0;
				ws.w_fftw = 0;
				ws.h_dst = 0;
				ws.w_dst = 0;
				ws.IsPeriodicPadding2Src = false;
		//		ws.IsCircularlyWrapForValidSame = false;
			}
			else{
				ws.h_fftw = find_closest_factor(h_src, ws.FFTW_FACTORS);
				ws.w_fftw = find_closest_factor(w_src, ws.FFTW_FACTORS);
				ws.h_dst = h_src - h_kernel + 1;
				ws.w_dst = w_src - w_kernel + 1;
				ws.IsPeriodicPadding2Src = false;
		//		ws.IsCircularlyWrapForValidSame = false;
			}
			break;
		
		case LINEAR_CONVOLUTION_VALID_UNPADDED: // NON PADDING
		case LINEAR_CORRELATION_VALID_UNPADDED: // NON PADDING
			if (ws.h_kernel > ws.h_src || ws.w_kernel > ws.w_src){
				//printf("Warning : The 'valid' convolution results in an empty matrix\n");
				ws.h_fftw = 0;
				ws.w_fftw = 0;
				ws.h_dst = 0;
				ws.w_dst = 0;
				ws.IsPeriodicPadding2Src = false;
		//		ws.IsCircularlyWrapForValidSame = false;
			}
			else{
				ws.h_fftw = h_src;
				ws.w_fftw = w_src;
				ws.h_dst = h_src - h_kernel + 1;
				ws.w_dst = w_src - w_kernel + 1;
				ws.IsPeriodicPadding2Src = false;
		//		ws.IsCircularlyWrapForValidSame = false;
			}
			break;

			// *********************************************
			// Valid&Same Linear convolution or correlation
			// *********************************************
		case LINEAR_CONVOLUTION_VALIDSAME: // PADDING
		case LINEAR_CORRELATION_VALIDSAME: // PADDING
			if (ws.h_kernel > ws.h_src || ws.w_kernel > ws.w_src){
				//printf("Warning : The 'validSame' convolution results in an empty matrix\n");
				ws.h_fftw = 0;
				ws.w_fftw = 0;
				ws.h_dst = 0;
				ws.w_dst = 0;
				ws.IsPeriodicPadding2Src = false;
			//	ws.IsCircularlyWrapForValidSame = false;
			}
			else{
				ws.h_dst = h_src;
				ws.w_dst = w_src;
				ws.h_fftw = find_closest_factor(h_src + h_kernel - 1, ws.FFTW_FACTORS);
				ws.w_fftw = find_closest_factor(w_src + w_kernel - 1, ws.FFTW_FACTORS);
				ws.IsPeriodicPadding2Src = true;
			//	ws.IsCircularlyWrapForValidSame = false;
			}
			break;

			// with circularly wrap to convolution or correlation result
		case LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP:
		case LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP:
			if (ws.h_kernel > ws.h_src || ws.w_kernel > ws.w_src){
				//printf("Warning : The 'validSame' convolution results in an empty matrix\n");
				ws.h_fftw = 0;
				ws.w_fftw = 0;
				ws.h_dst = 0;
				ws.w_dst = 0;
				ws.IsPeriodicPadding2Src = false;
			//	ws.IsCircularlyWrapForValidSame = true;
			}
			else{
				ws.h_dst = h_src;
				ws.w_dst = w_src;
				ws.h_fftw = find_closest_factor(h_src + h_kernel - 1, ws.FFTW_FACTORS);
				ws.w_fftw = find_closest_factor(w_src + w_kernel - 1, ws.FFTW_FACTORS);
				ws.IsPeriodicPadding2Src = false; // due to we will do suitable index control of the convolution/correlation result;
			//	ws.IsCircularlyWrapForValidSame = true;
			}
			break;

		case LINEAR_CONVOLUTION_VALIDSAME_UNPADDED: // NON PADDING
		case LINEAR_CORRELATION_VALIDSAME_UNPADDED: // NON PADDING
			if (ws.h_kernel > ws.h_src || ws.w_kernel > ws.w_src){
				//printf("Warning : The 'valid' convolution results in an empty matrix\n");
				ws.h_fftw = 0;
				ws.w_fftw = 0;
				ws.h_dst = 0;
				ws.w_dst = 0;
				ws.IsPeriodicPadding2Src = false;
	    //    	ws.IsCircularlyWrapForValidSame = false;
			}
			else{
				ws.h_fftw = h_src + h_kernel - 1;
				ws.w_fftw = w_src + w_kernel - 1;
				ws.h_dst = h_src;
				ws.w_dst = w_src;
				ws.IsPeriodicPadding2Src = true;
			//	ws.IsCircularlyWrapForValidSame = false;
			}
			break;


			// *********************************************
			// Full Circular convolution
			// *********************************************
		case CIRCULAR_CONVOLUTION_FULL: // PADDING
			// We here want to compute a circular convolution modulo h_dst, w_dst
			// These two variables must have been set before calling init_workscape !!
			ws.h_dst = h_src + h_kernel - 1;
			ws.w_dst = w_src + w_kernel - 1;
			ws.h_fftw = find_closest_factor(h_src + h_kernel - 1, ws.FFTW_FACTORS);
			ws.w_fftw = find_closest_factor(w_src + w_kernel - 1, ws.FFTW_FACTORS);
			ws.IsPeriodicPadding2Src = false;
		//	ws.IsCircularlyWrapForValidSame = false;
			break;

		case CIRCULAR_CONVOLUTION_FULL_UNPADDED: // NO PADDING
			// We here want to compute a circular convolution modulo h_dst, w_dst
			// These two variables must have been set before calling init_workscape !!
			ws.h_dst = h_src + h_kernel - 1;
			ws.w_dst = w_src + w_kernel - 1;
			ws.h_fftw = h_src + h_kernel - 1;
			ws.w_fftw = w_src + w_kernel - 1;
			ws.IsPeriodicPadding2Src = false;
		//	ws.IsCircularlyWrapForValidSame = false;
			break;

			// *********************************************
			// Same Circular convolution
			// *********************************************
		case CIRCULAR_CONVOLUTION_SAME: // PADDING
			// Circular convolution with optimal sizes
			ws.h_dst = std::max<int>(h_src, h_kernel);
			ws.w_dst = std::max<int>(w_src, w_kernel);
			ws.h_fftw = find_closest_factor(ws.h_dst, ws.FFTW_FACTORS);
			ws.w_fftw = find_closest_factor(ws.w_dst, ws.FFTW_FACTORS);
			ws.IsPeriodicPadding2Src = false;
		//	ws.IsCircularlyWrapForValidSame = false;
			break;

		case CIRCULAR_CONVOLUTION_SAME_UNPADDED: // NON PADDING
			// Circular convolution with optimal sizes
			ws.h_dst = std::max<int>(h_src, h_kernel);
			ws.w_dst = std::max<int>(w_src, w_kernel);
			ws.h_fftw = ws.h_dst;
			ws.w_fftw = ws.w_dst;
			ws.IsPeriodicPadding2Src = false;
		//	ws.IsCircularlyWrapForValidSame = false;
			break;


		default:
			ws.IsPeriodicPadding2Src = false;
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_SAME_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALID_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME_UNPADDED \n");

			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_FULL_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_SAME_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALID_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_UNPADDED \n");

			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME_UNPADDED \n");

			printf("   - LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP \n");

		} // end of switch-case

		int Num_FFT_Double_Inputs = ws.h_fftw * ws.w_fftw;
		// preserve a block of memory for array of src
		ws.in_src = new double[Num_FFT_Double_Inputs];
		// the non-redundant complex outputs
		// where the division is rounded down
		int Num_FFT_Complex_Outputs = ws.h_fftw * (ws.w_fftw / 2 + 1); 
		ws.out_src = (double*)fftw_malloc(sizeof(fftw_complex) * Num_FFT_Complex_Outputs);

		// preserve a block of memory for array of kernel
		ws.in_kernel = new double[Num_FFT_Double_Inputs];
		ws.out_kernel = (double*)fftw_malloc(sizeof(fftw_complex) * Num_FFT_Complex_Outputs);

		// normalized IFFT, i.e., normalized raw convolution, without indices selection
		ws.dst_ifft = new double[Num_FFT_Double_Inputs];  
		// final correlation/convolution corresponding to CONVOLUTION MODE, i.e., the selected elements of the normalized "dst_ifft"
		ws.dst_convolution = new double[ws.h_dst * ws.w_dst]; 

		// Initialization of the plans
		// 2-D array in row-major order, dim 0th is in the direction of height, 
		// and dim 1st is in the direction of width, 
		ws.p_forw_src =    fftw_plan_dft_r2c_2d(ws.h_fftw, ws.w_fftw, ws.in_src, (fftw_complex*)ws.out_src, ws.fftw_flag);
		ws.p_forw_kernel = fftw_plan_dft_r2c_2d(ws.h_fftw, ws.w_fftw, ws.in_kernel, (fftw_complex*)ws.out_kernel, ws.fftw_flag);

		// The backward FFT takes ws.out_kernel as input !!
		ws.p_back = fftw_plan_dft_c2r_2d(ws.h_fftw, ws.w_fftw, (fftw_complex*)ws.out_kernel, ws.dst_ifft, ws.fftw_flag);
	}

	void update_workspace_few(FFT_Workspace & ws){
		switch (ws.mode){
		case 0:
			ws.s_mode = "LINEAR_CONVOLUTION_FULL"; break;
		case 10:
			ws.s_mode = "LINEAR_CONVOLUTION_FULL_UNPADDED"; break;
		case 1:
			ws.s_mode = "LINEAR_CONVOLUTION_SAME";
			break;
		case 11:
			ws.s_mode = "LINEAR_CONVOLUTION_SAME_UNPADDED";
			break;
		case 2:
			ws.s_mode = "LINEAR_CONVOLUTION_VALID";
			break;
		case 12:
			ws.s_mode = "LINEAR_CONVOLUTION_VALID_UNPADDED";
			break;
		case 3:
			ws.s_mode = "LINEAR_CONVOLUTION_VALIDSAME";
			break;
		case 13:
			ws.s_mode = "LINEAR_CONVOLUTION_VALIDSAME_UNPADDED";
			break;
		case 4:
			ws.s_mode = "LINEAR_CORRELATION_FULL";
			break;
		case 14:
			ws.s_mode = "LINEAR_CORRELATION_FULL_UNPADDED";
			break;
		case 5:
			ws.s_mode = "LINEAR_CORRELATION_SAME";
			break;
		case 15:
			ws.s_mode = "LINEAR_CORRELATION_SAME_UNPADDED";
			break;
		case 6:
			ws.s_mode = "LINEAR_CORRELATION_VALID";
			break;
		case 16:
			ws.s_mode = "LINEAR_CORRELATION_VALID_UNPADDED";
			break;
		case 7:
			ws.s_mode = "LINEAR_CORRELATION_VALIDSAME";
			break;
		case 17:
			ws.s_mode = "LINEAR_CORRELATION_VALIDSAME_UNPADDED";
			break;
		case 8:
			ws.s_mode = "CIRCULAR_CONVOLUTION_FULL";
			break;
		case 18:
			ws.s_mode = "CIRCULAR_CONVOLUTION_FULL_UNPADDED";
			break;
		case 9:
			ws.s_mode = "CIRCULAR_CONVOLUTION_SAME";
			break;
		case 19:
			ws.s_mode = "CIRCULAR_CONVOLUTION_SAME_UNPADDED";
			break;
		case 20:
			ws.s_mode = "LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP";
			break;
		case 21:
			ws.s_mode = "LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP";
			break;
		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_SAME_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALID_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME_UNPADDED \n");

			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_FULL_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_SAME_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALID_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_UNPADDED \n");

			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME_UNPADDED \n");

			printf("   - LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP \n");
		}


		switch (ws.mode){
			// Here if no explicit "UNPADDED" occurs, we always mean PADDING as default.
			// i.e., here LINEAR_FULL == LINEAR_FULL_PADDED

		case LINEAR_CONVOLUTION_FULL: // PADDING
		case LINEAR_CORRELATION_FULL: // PADDING
		case LINEAR_CONVOLUTION_FULL_UNPADDED: // NON PADDING
		case LINEAR_CORRELATION_FULL_UNPADDED: // NON PADDING
		case LINEAR_CONVOLUTION_SAME_UNPADDED: // NO PADDING
		case LINEAR_CORRELATION_SAME_UNPADDED: // NO PADDING
		case LINEAR_CONVOLUTION_SAME: // PADDING
		case LINEAR_CORRELATION_SAME: // PADDING
		case LINEAR_CONVOLUTION_VALID: // PADDING
		case LINEAR_CORRELATION_VALID: // PADDING
		case LINEAR_CONVOLUTION_VALID_UNPADDED: // NON PADDING
		case LINEAR_CORRELATION_VALID_UNPADDED: // NON PADDING
		case LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP: // NON PADDING
		case LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP: // NON PADDING
		case CIRCULAR_CONVOLUTION_FULL: // PADDING
		case CIRCULAR_CONVOLUTION_FULL_UNPADDED: // NO PADDING
		case CIRCULAR_CONVOLUTION_SAME: // PADDING
		case CIRCULAR_CONVOLUTION_SAME_UNPADDED: // NON PADDING
			ws.IsPeriodicPadding2Src = false;
			//	ws.IsCircularlyWrapForValidSame = false;
			break;

			
			// *********************************************
			// Valid&Same Linear convolution or correlation
			// *********************************************
		case LINEAR_CONVOLUTION_VALIDSAME: // PADDING
		case LINEAR_CORRELATION_VALIDSAME: // PADDING
		case LINEAR_CONVOLUTION_VALIDSAME_UNPADDED: // NON PADDING
		case LINEAR_CORRELATION_VALIDSAME_UNPADDED: // NON PADDING
			if (ws.h_kernel > ws.h_src || ws.w_kernel > ws.w_src){
				//printf("Warning : The 'validSame' convolution results in an empty matrix\n");
				ws.IsPeriodicPadding2Src = false;
				//	ws.IsCircularlyWrapForValidSame = false;
			}
			else{
				ws.IsPeriodicPadding2Src = true;
				//	ws.IsCircularlyWrapForValidSame = false;
			}
			break;

		default:
			ws.IsPeriodicPadding2Src = false;
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_SAME_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALID_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME_UNPADDED \n");

			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_FULL_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_SAME_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALID_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_UNPADDED \n");

			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME_UNPADDED \n");

			printf("   - LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP \n");

		} // end of switch-case
	}

	void clear_workspace(FFT_Workspace & ws)
	{
		delete[] ws.FFTW_FACTORS;
		delete[] ws.in_src;
		fftw_free((fftw_complex*)ws.out_src);
		delete[] ws.in_kernel;
		fftw_free((fftw_complex*)ws.out_kernel);

		delete[] ws.dst_ifft;
		delete[] ws.dst_convolution;

		// Destroy the plans
		fftw_destroy_plan(ws.p_forw_src);
		fftw_destroy_plan(ws.p_forw_kernel);
		fftw_destroy_plan(ws.p_back);
	}

	// Only do FFT to src, instead of convolution or correlation
	// 指针常量－－指针本身是常量，指向的地址不可以变化,但是指向的地址所对应的内容可以变化
	// 常量指针－－e.g., 指向字符串常量，所指向的字符串内容不能变，但是指向的地址可以变化
	void fftw_src_fft(FFT_Workspace &ws, double * src){
		clock_t t = clock();

		double * ptr = 0, *ptr_end = 0;
		int Num_FFT_Double_Inputs = ws.h_fftw * ws.w_fftw;

		// Reset the content of ws.in_src as all zeros, for zero-padding of FFT
		for (ptr = ws.in_src, ptr_end = ws.in_src + Num_FFT_Double_Inputs; ptr != ptr_end; ++ptr)
			*ptr = 0.0;

		// two kinds of padding to the input signal: 1) just zero-padding ; 2) padding via circular shifting
		if (ws.IsPeriodicPadding2Src){ // In the case of ValidSame, we do padding to input signal via circular shifting,
			// making the padded input signal have the length of legth_for_fft_ifft
			for (int i = 0; i < ws.h_fftw; ++i)
				for (int j = 0; j < ws.w_fftw; ++j)
					ws.in_src[i*ws.w_fftw + j] = src[(i%ws.h_src)*ws.w_src + (j%ws.w_src)];
		}
		else{ // just do zero-padding to the input signal
			for (int i = 0; i < ws.h_src; ++i)
				for (int j = 0; j < ws.w_src; ++j)
					ws.in_src[i*ws.w_fftw + j] = src[i*ws.w_src + j];
		}

		// And we compute the packed FFT to src
		fftw_execute(ws.p_forw_src);

		if (ws.showSeconds){
			t = clock() - t;
			std::cout << "Function fftw_src_fft() in mode of " << ws.s_mode;
			printf(" took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
		}
	}


	 // only do FFT to Kernel, instead of convolution or correlation
	// For convolution we DO NOT flip the kernel;
	// But for correlation we HAVE TO flip the kernel
	// Note: So the reverse thing has been considered in this function,
	// and all others should not involve the flipping any more.
	// Otherwise, double flipping equals NO-Flipping, thus is wrong.
	void fftw_kernel_fft(FFT_Workspace &ws, double * kernel){
		clock_t t = clock();

		// to decide whether to flip the kernel or not
		switch (ws.mode){
		// *******************************************************
		// for linear convolution, Do NOt flip kernel
		// *******************************************************
		    // * full
		case	LINEAR_CONVOLUTION_FULL:
		case	LINEAR_CONVOLUTION_FULL_UNPADDED:
			// * same
		case	LINEAR_CONVOLUTION_SAME:
		case	LINEAR_CONVOLUTION_SAME_UNPADDED:
			// * valid
		case	LINEAR_CONVOLUTION_VALID:
		case	LINEAR_CONVOLUTION_VALID_UNPADDED:
			// * valid&same
		case	LINEAR_CONVOLUTION_VALIDSAME:
		case    LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP:
		case	LINEAR_CONVOLUTION_VALIDSAME_UNPADDED:
		
			if (ws.showSeconds)
				std::cout << "For linear convolution, DO NOT reverse kernel" << std::endl;
			break;


		// *******************************************************
		// for linear correlation, Flip kernel
		// *******************************************************
			// * full
		case LINEAR_CORRELATION_FULL:
		case LINEAR_CORRELATION_FULL_UNPADDED:
			// * same
		case LINEAR_CORRELATION_SAME:
		case LINEAR_CORRELATION_SAME_UNPADDED:
			// * valid
		case LINEAR_CORRELATION_VALID:
		case LINEAR_CORRELATION_VALID_UNPADDED:
			// * valid&same
		case LINEAR_CORRELATION_VALIDSAME:
		case LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP:
		case LINEAR_CORRELATION_VALIDSAME_UNPADDED:
			if (ws.showSeconds)
				std::cout << "For linear correlation, DO reverse kernel" << std::endl;
			std::reverse(kernel, kernel + ws.h_kernel*ws.w_kernel);
			break;

			// *******************************************************
			// for circular convolution, Do NOt flip kernel
			// *******************************************************
		case CIRCULAR_CONVOLUTION_SAME:
		case CIRCULAR_CONVOLUTION_FULL:
		case CIRCULAR_CONVOLUTION_SAME_UNPADDED:
		case CIRCULAR_CONVOLUTION_FULL_UNPADDED:
			if (ws.showSeconds)
				std::cout << "For circular convolution, DO NOT reverse kernel" << std::endl;
			break;

		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_SAME_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALID_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME_UNPADDED \n");

			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_FULL_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_SAME_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALID_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_UNPADDED \n");

			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME_UNPADDED \n");

			printf("   - LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP \n");
	} // end of switch-case 


		double * ptr = 0, *ptr_end = 0;
		int Num_FFT_Double_Inputs = ws.h_fftw * ws.w_fftw;

		// Reset the content of ws.in_kernel as all zeros, for zero-padding of FFT
		for (ptr = ws.in_kernel, ptr_end = ws.in_kernel + Num_FFT_Double_Inputs; ptr != ptr_end; ++ptr)
			*ptr = 0.0;

		// Only one kind of padding to the kernel, i.e., zero-padding
		// Then we do zero-padding to the kernel signal
		for (int i = 0; i < ws.h_kernel; ++i)
			for (int j = 0; j < ws.w_kernel; ++j)
				ws.in_kernel[i*ws.w_fftw + j] += kernel[i*ws.w_kernel + j];

		// And we compute the packed FFT to kernel
		fftw_execute(ws.p_forw_kernel);

		if (ws.showSeconds){
			t = clock() - t;
			std::cout << "Function fftw_kernel_fft() in mode of " << ws.s_mode;
			printf(" took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
		}
	}

	// to calculate convolution or correlation from two FFT signal, not the original time-domain or space domain signals
	void fftw_circular_convolution_from_FFT_Signal(FFT_Workspace &ws){
		clock_t t = clock();
		double * ptr = 0, *ptr_end = 0, *ptr2 = 0;
		int Num_FFT_Double_Inputs = ws.h_fftw * ws.w_fftw;

		// Compute the element-wise product on the packed terms
		// Let's put the element wise products to ws.in_kernel
		int Num_FFT_Complex_Outputs = ws.h_fftw * (ws.w_fftw / 2 + 1);
		double re_src, im_src,  // real part
			re_kernel, im_kernel; // imaginary part
		
		for (ptr = ws.out_src, ptr2 = ws.out_kernel, ptr_end = ws.out_src + 2 * Num_FFT_Complex_Outputs; ptr != ptr_end; ++ptr, ++ptr2)
		{
			re_src = *ptr;
			im_src = *(++ptr);
			re_kernel = *ptr2;
			im_kernel = *(++ptr2);
			*(ptr2 - 1) = re_src * re_kernel - im_src * im_kernel; // real part
			*ptr2 = re_src * im_kernel + im_src * re_kernel; // imaginary part
		}

		// Compute the backward FFT
		// Careful, The backward FFT does not preserve the output
		fftw_execute(ws.p_back);
		// Scale the transform
		for (ptr = ws.dst_ifft, ptr_end = ws.dst_ifft + Num_FFT_Double_Inputs; ptr != ptr_end; ++ptr)
			*ptr /= double(Num_FFT_Double_Inputs);

		if (ws.showSeconds){
			t = clock() - t;
			std::cout << "Function fftw_circular_convolution_from_FFT_Signal() in mode of " << ws.s_mode;
			printf(" took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
		}
}


	// Compute the circular convolution of src and kernel modulo ws.h_fftw, ws.w_fftw
	// using the Fast Fourier Transform (FFT)
	// The result is in ws.dst
	void fftw_circular_convolution(FFT_Workspace &ws, double * src, double * kernel)
	{
		clock_t t = clock();
	
		double * ptr = 0, *ptr_end = 0, *ptr2 = 0;
		// do FFT to src
		fftw_src_fft(ws, src);
		// do FFT to  kernel
		fftw_kernel_fft(ws, kernel);

		fftw_circular_convolution_from_FFT_Signal(ws);

		/*
		int Num_FFT_Double_Inputs = ws.h_fftw * ws.w_fftw;


		// Compute the element-wise product on the packed terms
		// Let's put the element wise products to ws.in_kernel
		int Num_FFT_Complex_Outputs = ws.h_fftw * (ws.w_fftw / 2 + 1);
		double re_src, im_src,  // real part
			   re_kernel, im_kernel; // imaginary part
		for (ptr = ws.out_src, ptr2 = ws.out_kernel, ptr_end = ws.out_src + 2 * Num_FFT_Complex_Outputs; ptr != ptr_end; ++ptr, ++ptr2)
		{
			re_src = *ptr;
			im_src = *(++ptr);
			re_kernel = *ptr2;
			im_kernel = *(++ptr2);
			*(ptr2 - 1) = re_src * re_kernel - im_src * im_kernel; // real part
			*ptr2 = re_src * im_kernel + im_src * re_kernel; // imaginary part
		}

		// Compute the backward FFT
		// Careful, The backward FFT does not preserve the output
		fftw_execute(ws.p_back);
		// Scale the transform
		for (ptr = ws.dst_ifft, ptr_end = ws.dst_ifft + Num_FFT_Double_Inputs; ptr != ptr_end; ++ptr)
			*ptr /= double(Num_FFT_Double_Inputs);

		if (ws.showSeconds){
			t = clock() - t;
			std::cout << "Function fftw_circular_convolution() in mode of " << ws.s_mode;
			printf(" took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
		}
		*/

		// That's it !
	}


	// /* function memcpy (#include <cstring>)
	//  * Syntax: void * memcpy ( void * destination, const void * source, size_t num ); // Copy block of memory
    //  * Copies the values of num bytes from the location pointed by source directly to the memory block pointed by destination.
    //  * The underlying type of the objects pointed by both the source and destination pointers are irrelevant for this function; 
	//  * The result is a binary copy of the data.
	//  * The function does not check for any terminating null character in source - it always copies exactly num bytes.
	//  * To avoid overflows, the size of the arrays pointed by both the destination and source parameters, shall be at least num bytes, and should not overlap (for overlapping memory blocks, memmove is a safer approach).
	//  */


	// Before calling this function, we must have gotten the FFT signals of src and kernel,
    // i.e., the FFT signal of src has been saved into ws.out_src, 
	// and the FFT signal of kernel has been saved into ws.out_kernel;
	// Otherwise, we will encounter errors !!!
	void convolve_from_FFT_Signal(FFT_Workspace &ws){
		clock_t t = clock();
		
		if (ws.h_fftw <= 0 || ws.w_fftw <= 0)
			return;

		// Depending on the type of convolution one is looking for, we extract the appropriate part of the result from ws.dst_ifft
		int h_offset, w_offset;
		switch (ws.mode){

			// *******************************************************
			// for linear convolution, no need of flipping of kernel
			// *******************************************************
			// * full
		case	LINEAR_CONVOLUTION_FULL:
		case	LINEAR_CONVOLUTION_FULL_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Full Linear convolution
			// Here we just keep the first [0:h_dst-1 ; 0:w_dst-1] elements of ws.dst_ifft (NOTE: dst_ifft only has real numbers )
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[i*ws.w_fftw], ws.w_dst*sizeof(double));
			break;

			// * same
		case	LINEAR_CONVOLUTION_SAME:
		case	LINEAR_CONVOLUTION_SAME_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Same linear convolution
			// Here we just keep the first [h_filter/2 : h_filter/2+h_dst-1 ; w_filter/2 : w_filter/2+w_dst-1] elements of ws.dst_ifft
			h_offset = int(ws.h_kernel / 2.0);
			w_offset = int(ws.w_kernel / 2.0);
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;

			// * valid
		case	LINEAR_CONVOLUTION_VALID:
		case	LINEAR_CONVOLUTION_VALID_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Valid linear convolution
			// Here we just take [h_dst x w_dst] elements starting at [h_kernel-1;w_kernel-1]
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;

			// * valid&same, without circularly wrap
		case	LINEAR_CONVOLUTION_VALIDSAME:
		case	LINEAR_CONVOLUTION_VALIDSAME_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Valid&Same linear convolution
			// Here we just take h_dst x w_dst elements starting at the point of (h_kernel-1, w_kernel-1);
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;

			// * valid&same, with circularly wrap
		case LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP: // circularly wrap

			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Valid&Same linear convolution
			// Here we just take h_dst x w_dst elements starting at the point of (h_kernel-1, w_kernel-1);
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;

			// circularly wrap, for two column blocks 
			// e.g., h_kernel = w_kernel = 8, we select the the first 7-column block and the last 7-column block, 
			// and send the sum of the two blocks to the first 7-column block.
			for (int r = 0; r < ws.h_fftw; ++ r)
				for (int c = 0; c < w_offset; ++ c)
					// ws.w_fftw*r + c : first 7-column block.
					// ws.w_fftw*r + (ws.w_dst + c) : last 7-column block.
					ws.dst_ifft[ws.w_fftw*r + c] += ws.dst_ifft[ws.w_fftw*r + (ws.w_dst + c)];

			// circularly wrap, for two row blocks 
			// e.g., e.g., h_kernel = w_kernel = 8, we select the the first 7-row block and the last 7-row block, 
			// and send the sum of the two row blocks to the first 7-row block.
			for (int r = 0; r < h_offset; ++r)
				for (int c = 0; c < ws.w_fftw; ++c)
					 //  ws.w_fftw*r + c : first 7-row block.
			         // ws.w_fftw*(r + ws.h_fftw) + c :the last 7-row block
					 ws.dst_ifft[ws.w_fftw*r + c] += ws.dst_ifft[ws.w_fftw*(r + ws.h_dst) + c];
			

			// do the subscripts control of the convolution results
			for (int r = 0; r < ws.h_dst; ++r)
				for (int c = 0; c < ws.w_dst; ++c)
					ws.dst_convolution[r* ws.w_dst + c] = ws.dst_ifft[r* ws.w_fftw + c];
			
			break;

			// *******************************************************
			// for linear correlation, with need of flipping of kernel
			// *******************************************************

			// * full
		case LINEAR_CORRELATION_FULL:
		case LINEAR_CORRELATION_FULL_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Full Linear correlation
			// Here we just keep the first [0:h_dst-1 ; 0:w_dst-1] elements of ws.dst_ifft (NOTE: dst_ifft only has real numbers )
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[i*ws.w_fftw], ws.w_dst*sizeof(double));
			break;

			// * same
		case LINEAR_CORRELATION_SAME:
		case LINEAR_CORRELATION_SAME_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Same linear correlation
			// Here we just keep the first [h_filter/2 : h_filter/2+h_dst-1 ; w_filter/2 : w_filter/2+w_dst-1] elements of ws.dst_ifft
			h_offset = int(ws.h_kernel / 2.0);
			w_offset = int(ws.w_kernel / 2.0);
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;

			// * valid
		case LINEAR_CORRELATION_VALID:
		case LINEAR_CORRELATION_VALID_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Valid linear correlation
			// Here we just take [h_dst x w_dst] elements starting at [h_kernel-1;w_kernel-1]
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;

			// * valid&same, without circularly wrap
		case LINEAR_CORRELATION_VALIDSAME:
		case LINEAR_CORRELATION_VALIDSAME_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Valid&Same linear correlation
			// Here we just take h_dst x w_dst elements starting at the point of (h_kernel-1, w_kernel-1);
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;


			// * valid&same, with circularly wrap
		case LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP: // circularly wrap
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// Valid&Same linear convolution
			// Here we just take h_dst x w_dst elements starting at the point of (h_kernel-1, w_kernel-1);
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;

			// circularly wrap, for two column blocks 
			// e.g., h_kernel = w_kernel = 8, we select the the first 7-column block and the last 7-column block, 
			// and send the sum of the two blocks to the first 7-column block.
			for (int r = 0; r < ws.h_fftw; ++r)
				for (int c = 0; c < w_offset; ++c)
					// ws.w_fftw*r + c : first 7-column block.
					// ws.w_fftw*r + (ws.w_dst + c) : last 7-column block.
					ws.dst_ifft[ws.w_fftw*r + c] += ws.dst_ifft[ws.w_fftw*r + (ws.w_dst + c)];

			// circularly wrap, for two row blocks 
			// e.g., e.g., h_kernel = w_kernel = 8, we select the the first 7-row block and the last 7-row block, 
			// and send the sum of the two row blocks to the first 7-row block.
			for (int r = 0; r < h_offset; ++r)
				for (int c = 0; c < ws.w_fftw; ++c)
					//  ws.w_fftw*r + c : first 7-row block.
					// ws.w_fftw*(r + ws.h_fftw) + c :the last 7-row block
					ws.dst_ifft[ws.w_fftw*r + c] += ws.dst_ifft[ws.w_fftw*(r + ws.h_dst) + c];


			// do the subscripts control of the convolution results
			for (int r = 0; r < ws.h_dst; ++r)
				for (int c = 0; c < ws.w_dst; ++c)
					ws.dst_convolution[r* ws.w_dst + c] = ws.dst_ifft[r* ws.w_fftw + c];
			break;

			// *******************************************************
			// for circular convolution
			// *******************************************************
		case CIRCULAR_CONVOLUTION_SAME:
		case CIRCULAR_CONVOLUTION_FULL:
		case CIRCULAR_CONVOLUTION_SAME_UNPADDED:
		case CIRCULAR_CONVOLUTION_FULL_UNPADDED:
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution_from_FFT_Signal(ws);
			// then, We copy the first [0:h_dst-1 ; 0:w_dst-1] elements of ws.dst_ifft for Circular convolution
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[i*ws.w_fftw], ws.w_dst*sizeof(double));
			break;

		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_SAME_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALID_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME_UNPADDED \n");

			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_FULL_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_SAME_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALID_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_UNPADDED \n");

			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME_UNPADDED \n");

			printf("   - LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP \n");
		} // end of switch-case 

		if (ws.showSeconds){
			t = clock() - t;
			std::cout << "Function convolve() in mode of " << ws.s_mode;
			printf(" took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
		}
		// that is it !!
	}

	void convolve(FFT_Workspace &ws, double * src, double * kernel){
		clock_t t = clock();

		if (ws.h_fftw <= 0 || ws.w_fftw <= 0)
			return;

		// Depending on the type of convolution one is looking for, we extract the appropriate part of the result from ws.dst_ifft
		int h_offset, w_offset;
		switch (ws.mode){

		// *******************************************************
		// for linear convolution, no need of flipping of kernel
		// *******************************************************
		// * full
		case	LINEAR_CONVOLUTION_FULL:
		case	LINEAR_CONVOLUTION_FULL_UNPADDED:
			// Compute the circular convolution
			fftw_circular_convolution(ws, src, kernel);
			// Full Linear convolution
			// Here we just keep the first [0:h_dst-1 ; 0:w_dst-1] elements of ws.dst_ifft (NOTE: dst_ifft only has real numbers )
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[i*ws.w_fftw], ws.w_dst*sizeof(double));
			break;

		// * same
		case	LINEAR_CONVOLUTION_SAME:
		case	LINEAR_CONVOLUTION_SAME_UNPADDED:
			// Compute the circular convolution
			fftw_circular_convolution(ws, src, kernel);
			// Same linear convolution
			// Here we just keep the first [h_filter/2 : h_filter/2+h_dst-1 ; w_filter/2 : w_filter/2+w_dst-1] elements of ws.dst_ifft
			h_offset = int(ws.h_kernel / 2.0);
			w_offset = int(ws.w_kernel / 2.0);
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;
			
		// * valid
		case	LINEAR_CONVOLUTION_VALID:
		case	LINEAR_CONVOLUTION_VALID_UNPADDED:
			// Compute the circular convolution
			fftw_circular_convolution(ws, src, kernel);
			// Valid linear convolution
			// Here we just take [h_dst x w_dst] elements starting at [h_kernel-1;w_kernel-1]
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;
			
		// * valid&same
		case	LINEAR_CONVOLUTION_VALIDSAME:
		case	LINEAR_CONVOLUTION_VALIDSAME_UNPADDED:
			// Compute the circular convolution
			fftw_circular_convolution(ws, src, kernel);
			// Valid&Same linear convolution
			// Here we just take h_dst x w_dst elements starting at the point of (h_kernel-1, w_kernel-1);
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;


			// * valid&same, with circularly wrap
		case LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP: // circularly wrap
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution(ws, src, kernel);
			// Valid&Same linear convolution
			// Here we just take h_dst x w_dst elements starting at the point of (h_kernel-1, w_kernel-1);
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;

			// circularly wrap, for two column blocks 
			// e.g., h_kernel = w_kernel = 8, we select the the first 7-column block and the last 7-column block, 
			// and send the sum of the two blocks to the first 7-column block.
			for (int r = 0; r < ws.h_fftw; ++r)
				for (int c = 0; c < w_offset; ++c)
					// ws.w_fftw*r + c : first 7-column block.
					// ws.w_fftw*r + (ws.w_dst + c) : last 7-column block.
					ws.dst_ifft[ws.w_fftw*r + c] += ws.dst_ifft[ws.w_fftw*r + (ws.w_dst + c)];

			// circularly wrap, for two row blocks 
			// e.g., e.g., h_kernel = w_kernel = 8, we select the the first 7-row block and the last 7-row block, 
			// and send the sum of the two row blocks to the first 7-row block.
			for (int r = 0; r < h_offset; ++r)
				for (int c = 0; c < ws.w_fftw; ++c)
					//  ws.w_fftw*r + c : first 7-row block.
					// ws.w_fftw*(r + ws.h_fftw) + c :the last 7-row block
					ws.dst_ifft[ws.w_fftw*r + c] += ws.dst_ifft[ws.w_fftw*(r + ws.h_dst) + c];


			// do the subscripts control of the convolution results
			for (int r = 0; r < ws.h_dst; ++r)
				for (int c = 0; c < ws.w_dst; ++c)
					ws.dst_convolution[r* ws.w_dst + c] = ws.dst_ifft[r* ws.w_fftw + c];
			break;
		
			// *******************************************************
			// for linear correlation, with need of flipping of kernel
			// *******************************************************
        
		// * full
		case LINEAR_CORRELATION_FULL:
		case LINEAR_CORRELATION_FULL_UNPADDED:
			
			// compute the circular convolution
			   fftw_circular_convolution(ws, src, kernel);
			   // Note: the above function will call fuction fftw_kernel_fft()
			   //  and the reverse thing has been considered in function fftw_kernel_fft(),
			   // all others should not involve the flipping any more.
			   // Otherwise, double flipping equals NO-Flipping, thus is wrong.

			// Full Linear correlation
			// Here we just keep the first [0:h_dst-1 ; 0:w_dst-1] elements of ws.dst_ifft (NOTE: dst_ifft only has real numbers )
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[i*ws.w_fftw], ws.w_dst*sizeof(double));
			break;

		// * same
		case LINEAR_CORRELATION_SAME:
		case LINEAR_CORRELATION_SAME_UNPADDED:
			// compute the circular convolution
				fftw_circular_convolution(ws, src, kernel);
				// Note: the above function will call fuction fftw_kernel_fft()
				//  and the reverse thing has been considered in function fftw_kernel_fft(),
				// all others should not involve the flipping any more.
				// Otherwise, double flipping equals NO-Flipping, thus is wrong.
			
			// Same linear correlation
			// Here we just keep the first [h_filter/2 : h_filter/2+h_dst-1 ; w_filter/2 : w_filter/2+w_dst-1] elements of ws.dst_ifft
			h_offset = int(ws.h_kernel / 2.0);
			w_offset = int(ws.w_kernel / 2.0);
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;

		// * valid
		case LINEAR_CORRELATION_VALID:
		case LINEAR_CORRELATION_VALID_UNPADDED:
			// compute the circular convolution
				fftw_circular_convolution(ws, src, kernel);
				// Note: the above function will call fuction fftw_kernel_fft()
				//  and the reverse thing has been considered in function fftw_kernel_fft(),
				// all others should not involve the flipping any more.
				// Otherwise, double flipping equals NO-Flipping, thus is wrong.

			// Valid linear correlation
			// Here we just take [h_dst x w_dst] elements starting at [h_kernel-1;w_kernel-1]
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;
			
		// * valid&same
		case LINEAR_CORRELATION_VALIDSAME:
		case LINEAR_CORRELATION_VALIDSAME_UNPADDED:
             // compute the circular convolution
				fftw_circular_convolution(ws, src, kernel);
				// Note: the above function will call fuction fftw_kernel_fft()
				//  and the reverse thing has been considered in function fftw_kernel_fft(),
				// all others should not involve the flipping any more.
				// Otherwise, double flipping equals NO-Flipping, thus is wrong.

			// Valid&Same linear correlation
			// Here we just take h_dst x w_dst elements starting at the point of (h_kernel-1, w_kernel-1);
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[(i + h_offset)*ws.w_fftw + w_offset], ws.w_dst*sizeof(double));
			break;
			
			// * valid&same, with circularly wrap
		case LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP: // circularly wrap
			// Compute the circular convolution of 2 FFT signals
			// i.e., first do wise multiplication of two FFT signals, and then do IFFT of the multiplication
			fftw_circular_convolution(ws, src, kernel);
			// Valid&Same linear convolution
			// Here we just take h_dst x w_dst elements starting at the point of (h_kernel-1, w_kernel-1);
			h_offset = ws.h_kernel - 1;
			w_offset = ws.w_kernel - 1;

			// circularly wrap, for two column blocks 
			// e.g., h_kernel = w_kernel = 8, we select the the first 7-column block and the last 7-column block, 
			// and send the sum of the two blocks to the first 7-column block.
			for (int r = 0; r < ws.h_fftw; ++r)
				for (int c = 0; c < w_offset; ++c)
					// ws.w_fftw*r + c : first 7-column block.
					// ws.w_fftw*r + (ws.w_dst + c) : last 7-column block.
					ws.dst_ifft[ws.w_fftw*r + c] += ws.dst_ifft[ws.w_fftw*r + (ws.w_dst + c)];

			// circularly wrap, for two row blocks 
			// e.g., e.g., h_kernel = w_kernel = 8, we select the the first 7-row block and the last 7-row block, 
			// and send the sum of the two row blocks to the first 7-row block.
			for (int r = 0; r < h_offset; ++r)
				for (int c = 0; c < ws.w_fftw; ++c)
					//  ws.w_fftw*r + c : first 7-row block.
					// ws.w_fftw*(r + ws.h_fftw) + c :the last 7-row block
					ws.dst_ifft[ws.w_fftw*r + c] += ws.dst_ifft[ws.w_fftw*(r + ws.h_dst) + c];


			// do the subscripts control of the convolution results
			for (int r = 0; r < ws.h_dst; ++r)
				for (int c = 0; c < ws.w_dst; ++c)
					ws.dst_convolution[r* ws.w_dst + c] = ws.dst_ifft[r* ws.w_fftw + c];
			break;

			// *******************************************************
			// for circular convolution
			// *******************************************************
		case CIRCULAR_CONVOLUTION_SAME:
		case CIRCULAR_CONVOLUTION_FULL:
		case CIRCULAR_CONVOLUTION_SAME_UNPADDED:
		case CIRCULAR_CONVOLUTION_FULL_UNPADDED:
			// first, to compute the circular convolution
			fftw_circular_convolution(ws, src, kernel);
			
			// then, We copy the first [0:h_dst-1 ; 0:w_dst-1] elements of ws.dst_ifft for Circular convolution
			for (int i = 0; i < ws.h_dst; ++i)
				memcpy(&ws.dst_convolution[i*ws.w_dst], &ws.dst_ifft[i*ws.w_fftw], ws.w_dst*sizeof(double));
			break;

		default:
			printf("Unrecognized convolution mode, possible modes are :\n");
			printf("   - LINEAR_CONVOLUTION_FULL \n");
			printf("   - LINEAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_SAME \n");
			printf("   - LINEAR_CONVOLUTION_SAME_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALID \n");
			printf("   - LINEAR_CONVOLUTION_VALID_UNPADDED \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME \n");
			printf("   - LINEAR_CONVOLUTION_VALIDSAME_UNPADDED \n");

			printf("   - LINEAR_CORRELATION_FULL \n");
			printf("   - LINEAR_CORRELATION_FULL_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_SAME \n");
			printf("   - LINEAR_CORRELATION_SAME_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALID \n");
			printf("   - LINEAR_CORRELATION_VALID_UNPADDED \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_UNPADDED \n");

			printf("   - CIRCULAR_CONVOLUTION_FULL \n");
			printf("   - CIRCULAR_CONVOLUTION_FULL_UNPADDED \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME \n");
			printf("   - CIRCULAR_CONVOLUTION_SAME_UNPADDED \n");

			printf("   - LINEAR_CONVOLUTION_VALIDSAME_CIRCULAR_WRAP \n");
			printf("   - LINEAR_CORRELATION_VALIDSAME_CIRCULAR_WRAP \n");
		} // end of switch-case 

		if (ws.showSeconds){
		t = clock() - t;
		std::cout << "Function convolve() in mode of " << ws.s_mode ;
		printf(" took %f seconds.\n", ((float)t) / CLOCKS_PER_SEC);
	}
}

	void displayConvolve(FFT_Workspace & ws){
		std::cout << "Convolution/Correlation mode is : " << ws.s_mode << std::endl;
		std::cout << "Convolution/Correlation 2-D result : " << std::endl;
		for (int i = 0; i < ws.h_dst; ++i){
			for (int j = 0; j < ws.w_dst; ++j){
				std::cout << ws.dst_convolution[i*ws.w_dst + j] << ", ";
			}
			std::cout << std::endl;
		}
	}

	void saveConvolve(FFT_Workspace & ws, std::string & filename){
		std::ofstream fout(filename,std::ios::app);
		if (!fout){
			std::cout << "File Not Opened" << std::endl;
		}
		fout << "Input Signal : " << ws.h_src << " by " << ws.w_src << std::endl
		     << "FFT Convolution/Correlation mode is : " << ws.s_mode << std::endl
			 << "FFT Convolution/Correlation 2-D result : " << std::endl;
		for (int i = 0; i < ws.h_dst; ++i){
			for (int j = 0; j < ws.w_dst; ++j){
				fout << ws.dst_convolution[i*ws.w_dst + j] << ", ";
			}
			fout << std::endl;
		}

		fout << std::endl;
		fout.close();
	}

	void help(){
		std::cout << "There are two methods to calculate convolution or correlation.\n"
			<< "//*************************************************************************\n"
			<< "The first one is to directly call function convolve(ws, src, kernel).\n"
			<< "//*************************************************************************\n\n"
			<< "//*************************************************************************\n"
			<< "The second one is to continuously call 3 functions as below:\n"
			<< "fftw_src_fft(ws, src);\n"
			<< "fftw_kernel_fft(ws, kernel);\n"
			<< "convolve_from_FFT_Signal(ws);\n"
			<< "//*************************************************************************\n";
	}

} // end of namespace FFTW_Convolution