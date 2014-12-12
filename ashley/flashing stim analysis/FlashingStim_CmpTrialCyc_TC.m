mouse = '516';
date = '141003';
ImgFolder = '002+003+004';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
load('dataStructDFoverF.mat')
load('dataMasks.mat')

%% save timecourses for different masks along with relevant variables
siz = size(dataStructDFoverF.Cycles,2);

%dirTuning
for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    dirTuning_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.dirTuning);
    clear thisData
end

for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    FSfirstRsp_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.FSfirstRsp);
    clear thisData
end

for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    FScycMaxDF_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.FScycMaxDF{icyc});
    clear thisData
end

for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    FScycTargetRsp_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.FScycTargetRsp{icyc});
    clear thisData
end

dataTC.dirTuning = dirTuning_TC;
dataTC.FSfirstRsp = FSfirstRsp_TC;
dataTC.FScycMaxDF = FScycMaxDF_TC;
dataTC.FScycTargetRsp = FScycTargetRsp_TC;

save('dataTC.mat','dataTC')
