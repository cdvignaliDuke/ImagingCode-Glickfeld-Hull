iexp = 1;
naiveAudModDatasets_V1
rc = behavConstsNAM;

irun = 1;


%% analysis folder
try
    filedir = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,'two-photon imaging', expt(iexp).date, expt(iexp).runs(irun,:));
    cd(filedir);
catch
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging');
    cd(filedir)
    mkdir(date,ImgFolder)
    filedir = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(filedir);
end

fnout = filedir;

%% align data
data = double(data);

%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
clear data_down

% register
data_avg = mean(data_sub(:,:,101:200),3);
figure; imagesq(data_avg); colormap(gray)

    %save registration image
    save(fullfile(fnout,'regImg.mat'),'data_avg');

[out data_reg] = stackRegister(data_sub, data_avg);
clear data_sub

%% get dF/F 
%get trial F
oneS = double(input.frameRateHz);
trialStart = double(cell2mat(input.cLeverDown));
itiFrames = linspaceNDim(trialStart',trialStart'-(oneS-1),oneS);
nTrials = length(trialStart);

dataF = mean(mean(reshape(data_reg(:,:,itiFrames(:)),size(data_reg,1),size(data_reg,2),nTrials,oneS),4),3);
dataDF = bsxfun(@minus,data_reg,dataF);
dataDFoverF = bsxfun(@rdivide,dataDF,dataF);

DFoverFdown = stackGroupProject(dataDFoverF,3);
writetiff(DFoverFdown,fullfile(fnout,'DFoverFdown'))

%% get maxDF/F
maxDFoverF = max(dataDFoverF,[],3);
figure;imagesc(maxDFoverF);colormap gray;