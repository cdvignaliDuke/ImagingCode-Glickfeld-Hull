CD = ('Z:\2P imaging\Analysis\516\141003\002\FlashingStimAnalysis\Rsp2OriChange');
cd(CD);
data002 = readtiff('PrePostTargetSuccess.tif');
CD = ('Z:\2P imaging\Analysis\516\141003\003\FlashingStimAnalysis\Rsp2OriChange');
cd(CD);
data003 = readtiff('PrePostTargetSuccess.tif');
CD = ('Z:\2P imaging\Analysis\516\141003\004\FlashingStimAnalysis\Rsp2OriChange');
cd(CD);
data004 = readtiff('PrePostTargetSuccess.tif');

data = cat(3,data002,data003,data004);
clear data002 data003 data004;

SubNum = '516';
date = '141003';
mouse = '516';
ImgFolder = '002+003+004';
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2OriChange'];
cd(CD);

% writetiff(data,'PrePostTargetSuccess.tif');
data = readtiff('PrePostTargetSuccess.tif');

CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);
load('mworksAll.mat')

start = 1;
for i = 1:3
    siz = size(inputAll(i).trialOutcomeCell,2);
    TrialOutcome(1,start:start+siz-1) = inputAll(i).trialOutcomeCell;
    start = start+siz;
end
success_log = strcmp(TrialOutcome,'success');
success_ind = find(success_log == 1);
success_ind = success_ind(:,[1:38 40:end]);

tTrials = size(success_ind,2);

start = 1;
for i = 1:3
    siz = size(inputAll(i).tBlock2TrialNumber,2);
    Block2On(1,start:start+siz-1) = cell2mat(inputAll(i).tBlock2TrialNumber);
    start = start+siz;
end
Block2OnSuccess = Block2On(:,success_ind);

start = 1;
for i = 1:3
    siz = size(inputAll(i).tGratingDirectionDeg,2);
    targetDir(1,start:start+siz-1) = cell2mat(inputAll(i).tGratingDirectionDeg);
    start = start+siz;
end
targetDirSuccess = targetDir(:,success_ind);
Dir = unique(targetDirSuccess);
DirRound = floor(Dir);
L = 30;
%%
F_byTrial = zeros(264,1250,tTrials);
for itrial = 1:tTrials
    ind = (itrial-1).*L+1:itrial.*L;
    F_byTrial(:,:,itrial) = mean(data(:,:,ind),3);
end

F = mean(data,3);

dF_data = zeros(size(data));
for itrial = 1:tTrials
    start1 = L.*(itrial-1)+1;
    start2 = L.*itrial;
    dF_data(:,:,start1:start2) = bsxfun(@minus,data(:,:,start1:start2),F);
end

dFoverF_data = zeros(size(data));
for itrial = 1:tTrials
    start1 = L.*(itrial-1)+1;
    start2 = L.*itrial;
    dFoverF_data(:,:,start1:start2) = bsxfun(@rdivide,dF_data(:,:,start1:start2),F);
end

last_ind = zeros(1,15);
dF_TargetRspMax = zeros(size(dFoverF_data,1),size(dFoverF_data,2),tTrials);
for itrial = 1:tTrials
    last_ind = (L.*itrial)-14:L.*itrial;
    dF_TargetRspMax(:,:,itrial) = max(dFoverF_data(:,:,last_ind),[],3);
end
clear last_ind

data_avg = mean(dFoverF_data,3);
figure;imagesq(data_avg);colormap(gray)


bwout = imCellEditInteractive(data_avg);
mask_cellTargetAll = bwlabel(bwout);

data_TC = stackGetTimeCourses(dFoverF_data,mask_cellTargetAll);
figure; tcOffsetPlot(data_TC)