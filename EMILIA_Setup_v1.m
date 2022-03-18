close all; clear all; clc;                                                                            for Advanced_Setup = 1
% - EMILIA Advanced Setup [FIXED] ------------------------------------------------------------------------------------------   
    Decimate  = 'No';    %------------------------------------------------------  Decimate data                          ---
    Factor    = 2;       %------------------------------------------------------  Factor to decimate data                ---
    TimeAlign = 'Yes';   %------------------------------------------------------  Re-Align delayed outputs               ---
    Wavelet   = 'db45';  %------------------------------------------------------  Wavelet function                       --- 
    Distri    = 'Kernel';%------------------------------------------------------  Probability distribution               ---
    Dist_Res  = 1000;    %------------------------------------------------------  PDF frequency vector resolution        ---
    Max_Damp  = [];      %------------------------------------------------------  Test maximum damping (0.05 = 5%)       ---
    Min_Prob  = [];      %------------------------------------------------------  Test minimum probability threshold     ---
    Include_Z = 'No';    %------------------------------------------------------  Include vertical axis data to analysis ---
% --------------------------------------------------------------------------------------------------------------------------   
                                                                                                     end
% - EMILIA Analysis Setup --------------------------------------------------------------------------------------------------
    Data_File = '2D_Bridge_v2_ST00_26ch_100Hz.dat'; %---------------------------  Raw acceleration data file.            ---                    
    Geo_File  = '2D_Bridge_v2_Geo.xls';   %-------------------------------------  Geometry data.                         ---
    fs        = 100;     %------------------------------------------------------  Sampling frequency.                    ---
    Level     = 3;       %------------------------------------------------------  Decomposition level.                   ---
% --------------------------------------------------------------------------------------------------------------------------
                                                                                                                for Code = 1                                                                                                              
[Acc_Data,Displa_Data_RAW,iFr,iAm,fs,Level,Channels,Bayes_PD,PDF_Fun,f_vector,MEV_x,Test_Array] = EMILIA_core_v1...
(Data_File,Geo_File,fs,Level,Decimate,Factor,TimeAlign,Wavelet,Distri,Dist_Res,Max_Damp,Min_Prob,Include_Z);
                                                                                                                end
%% --- Hilbert Spectrum ----------------------------------------------------------------------------------------------------
    Ch_plot   = [];      %------------------------------------------------------  Channels to plot.                      ---
    Sq_plot   = 1:8;     %------------------------------------------------------  Subsequence to plot.                   ---
    Sp_skip   = 250;     %------------------------------------------------------  Samples to skip.                       ---
    freq_span = [0 10];  %------------------------------------------------------  Frequency span (Y axis).               ---
    time_span = [];      %------------------------------------------------------  Time span (X axis).                    ---
    Theme     = 'Blue';  %------------------------------------------------------  Background colour.                     ---
    Order     = 'Local'; %------------------------------------------------------  Sub sequences plot order.              ---
    Dist_Res  = 2;       %------------------------------------------------------  Scatter plot resolution.               ---    
    Mov_AVG   = 'Yes';   %------------------------------------------------------  Plot Moving averages.                  ---   
    Mov_WIN   = fs;      %------------------------------------------------------  Moving averages time window length.    ---  
    Mov_skip  = 25;      %------------------------------------------------------  Samples to skip for moving averages    ---
% -------------------------------------------------------------------------------------------------------------------------- 
                                                                                                                for Code = 1
EMILIA_Hilbert_Spectrum_v1(Acc_Data,iFr,iAm,fs,Level,Channels,Ch_plot,Sq_plot,Sp_skip,...
freq_span,time_span,Theme,Order,Dist_Res,Wavelet,Mov_AVG,Mov_WIN,Mov_skip);
                                                                                                                end
%% --- Probability Spectrum ------------------------------------------------------------------------------------------------
    Ch_plot   = [];      %------------------------------------------------------  Channels to plot.                      ---
    Sq_plot   = 1:8;     %------------------------------------------------------  Sub sequences to plot.                 ---
    Sp_skip   = 250;     %------------------------------------------------------  Samples to skip.                       ---
    Y_Axis    = 'Lin';   %------------------------------------------------------  Y axis mode ('Lin'/'Log').             ---
    Plot_KePD = 'Yes';   %------------------------------------------------------  Plot Kernel PDFs ('Yes'/'No').         ---
    Plot_BaPD = 'Yes';   %------------------------------------------------------  Plot Bayes Likehood PDFs ('Yes'/'No'). ---
    Plot_Test = 'Yes';   %------------------------------------------------------  Plot Test results (Yes/No).            ---
    fn_span   = [0 10];  % -----------------------------------------------------  Frequency span to plot (X axis).       ---
% --------------------------------------------------------------------------------------------------------------------------
                                                                                                                for Code = 1
EMILIA_Probability_Spectrum_v1(Ch_plot,Sq_plot,Sp_skip,iFr,iAm,Bayes_PD,...
PDF_Fun,Wavelet,Level,fs,Channels,f_vector,MEV_x,Test_Array,Min_Prob,...
Y_Axis,Plot_KePD,Plot_BaPD,Plot_Test,fn_span);
                                                                                                                end
%% --- 2D Mode Shape Animation ---------------------------------------------------------------------------------------------
    EMILIA_Seq  = 1;     %------------------------------------------------------   Sub sequence to animate.              ---
    IniSample   = 1000;  %------------------------------------------------------   Animation initial sample.             ---
    EndSample   = 2000;  %------------------------------------------------------   Animation end sample.                 ---
    Make_Video  = 'No';  %------------------------------------------------------   Option to make video.                 ---     
    Video_Name  = 'video_name.avi'; %-------------------------------------------   Video file type and name.             ---
% --------------------------------------------------------------------------------------------------------------------------
                                                                                                                for Code = 1
    EMILIA_2D_Shape_Animation_v1(Displa_Data_RAW,MEV_x,EMILIA_Seq,IniSample,EndSample,Geo_File,Make_Video,Video_Name);  
                                                                                                                end
%% --- 3D Mode Shape Animation ----------------------------------------------------------------------------------------------
    EMILIA_Seq  = 2;      %-----------------------------------------------------   Sub sequence to animate.              ---
    IniSample   = 1000;   %-----------------------------------------------------   Animation initial sample.             ---
    EndSample   = 2000;   %-----------------------------------------------------   Animation end sample.                 ---  
    Make_Video  = 'No';   %-----------------------------------------------------   Option to make video.                 ---     
    Video_Name  = 'video_name.avi'; %-------------------------------------------   Video file type and name.             ---
% --------------------------------------------------------------------------------------------------------------------------
                                                                                                                for Code = 1
    EMILIA_3D_Shape_Animation_v1(Displa_Data_RAW,MEV_x,EMILIA_Seq,IniSample,EndSample,Geo_File,Make_Video,Video_Name);
                                                                                                                end