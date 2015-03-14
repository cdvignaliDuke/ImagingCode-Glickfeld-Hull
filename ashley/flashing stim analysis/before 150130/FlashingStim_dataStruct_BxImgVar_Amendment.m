%% Amendment to dataStructVar for FlashingStim_2P_Frames ON&OFF times
date = '141209';
mouse = 'AW07';
SubNum = '607';
ImgFolder = '002';
time = '1744';

    % MWorks file
CD = ['Z:\data\' mouse '\MWorks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

    % save in
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);

%% add variables

dataStructVar.ONfr = input.nFramesOn;
dataStructVar.OFFfr = input.nFramesOff;

save('dataStructVar.mat','dataStructVar');