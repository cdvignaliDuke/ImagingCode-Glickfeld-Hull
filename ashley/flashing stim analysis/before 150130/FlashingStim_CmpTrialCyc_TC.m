mouse = 'AW07';
SubNum = '607';
date = '150121';
ImgFolder = '003';
CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
%run FlahsingStim_dataStruct.m to organize data first
load('dataMasks.mat')

%% save timecourses for different masks along with relevant variables
siz = size(dataStructDFoverF.Cycles,2);

%dirTuning
for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    dirTuning_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.dirTuning);
    clear thisData
end

%retTuning
for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    retTuning_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.retTuning);
    clear thisData
end


% others

% for icyc = 1:siz
%     thisData = dataStructDFoverF.cycData{icyc};
%     FSfirstRsp_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.FSfirstRsp);
%     clear thisData
% end

for icyc = 1:siz
    thisData = dataStructDFoverF.cycData{icyc};
    FScycCorrDF_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.FScycCorrDF{icyc});
    clear thisData
end

% for icyc = 1:siz
%     thisData = dataStructDFoverF.cycData{icyc};
%     FScycTargetRsp_TC{icyc} = stackGetTimeCourses(thisData,dataMasks.FScycTargetRsp{icyc});
%     clear thisData
% end

dataTC.dirTuning = dirTuning_TC;
dataTC.retTuning = retTuning_TC;
% dataTC.FSfirstRsp = FSfirstRsp_TC;
dataTC.FScycCorrDF = FScycCorrDF_TC;
% dataTC.FScycTargetRsp = FScycTargetRsp_TC;

save('dataTC.mat','dataTC')
