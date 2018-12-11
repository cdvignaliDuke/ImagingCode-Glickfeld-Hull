%CRP_analysis_overview
%script which has the general outline for all current analyses performed on CRP data 
% note: does not include the code written purely for behavioral analysis during training. 

%% initial processing
% PCA/ICA
% bx analysis
% TC extraction
main_CRP2;

%% isolate Ca Events and generate PSTHs
%also analyzes Ca Event amplitude and waveform BEFORE removing outlier events
main2_CRP;

%% summary plots for grand means

%plots Ca event and licking PSTHs
summary_PSTHs_CRP_JH;
%plot df/f TC grand means
summary_TCs_CRP_JH;

