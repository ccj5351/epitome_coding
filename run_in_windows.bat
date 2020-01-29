::This experiment is used to compare different cased of parameter set-up.
::Dataset is Edwin-less-try, including 10 images.

::set appPath=C:/Users/ccj/Desktop/CPlusPlusProjects_CCJ/Epitome_6_5/Release/Epitome_6_5.exe
set appPath=C:/OpenCVProjects_CCJ/Epitome_FFT/x64/Release/Epitome_FFT.exe
set DatabaseDir="C:/EpitomeDataset/Stevens-sub-images/Edwin-less-try"
set EpitomeResultDir="C:/EpitomeDataset/epi-recon/epi-FFT-Edwin-less-try"
set ReconsCompresDir="C:/EpitomeDataset/compress-epi-recon/epi-FFT_compre-Edwin-less-try"
set display="0"
set whatkindImgs="Edwin-less-try"
set newPatchSpacingNum="1"
set newPatchSpacing="4"
set fileNameforRowCol_Idx_YML="baseline"
rem set fileNameforRowCol_Idx_YML="1-0-050-2-1"
set IsGrayImgs="0"
set save_img_num="10"
set nthreads="4"
set max_omp_for_idx=4
set randomshifting="0"
set epitomeWidth="256"
set patchSideLengh="8"
set patchSpacing="4"
set numIteration="10"
set nameDifference="recon-"
set imgEpi_EncodeType=".bmp"
set overheadNum="20"
set s_flag_epi_learn_recon="2"
set errorImgBit="8"
set jpeg_quality="90"
set s_standard_encodeType="JPEG2K"
set errHistThres="2"
set s_quantizeType="NOTHING"
set compres_flag="0"
set JPEG2K_EXE_Base_Path="C:/Users/ccj/Desktop/CPlusPlusProjects_CCJ/openjpeg-2.1.0-win32-x86/bin"
set compres_flag_2_epi_recon="0"
set basePath_Idx="C:/Users/ccj/imageDatabase/bmp-images-epi-result/epi-6_Try_epi-sub-Edwin-less-full/1-0-30-2-1"
set pauseTime="50"




  FOR  %%W IN (90) DO (start %appPath% %DatabaseDir% %EpitomeResultDir% %ReconsCompresDir% %display% %whatkindImgs% %newPatchSpacingNum% %newPatchSpacing%^
 %fileNameforRowCol_Idx_YML% %IsGrayImgs% %save_img_num% %nthreads% %max_omp_for_idx% %randomshifting% %epitomeWidth% %patchSideLengh% %patchSpacing% %numIteration% %nameDifference% %imgEpi_EncodeType% %overheadNum%^
 %s_flag_epi_learn_recon% %errorImgBit% %%W %s_standard_encodeType% %errHistThres% %s_quantizeType% %compres_flag% %JPEG2K_EXE_Base_Path%^
 %compres_flag_2_epi_recon% %basePath_Idx%)