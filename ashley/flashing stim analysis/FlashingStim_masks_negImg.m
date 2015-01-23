mouse = 'AW07';
SubNum = '607';
date = '150121';
ImgFolder = '004';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructVar.mat');
load('dataMasks.mat');
load('dataTC.mat');
%run FlahsingStim_dataStruct.m to organize data first

%%
mask = dataMasks.retTuning;
negativeMask = mask == 0;

dataMasks.retTuningNegative = negativeMask;

for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    retTuningNeg_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.retTuningNegative);
    clear thisData
end

dataTC.retTuningNegative = retTuningNeg_TC;
save('dataTC.mat','dataTC')
