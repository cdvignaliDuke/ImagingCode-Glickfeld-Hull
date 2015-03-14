edit Load_SBXdataset_fast.m

data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub;
clear data_sub
%register to averaged frames
data_avg = mean(data(:,:,5000:5090),3);
figure; imagesq(data_avg); colormap(gray)

[out data_reg] = stackRegister(data, data_avg);
clear data
data = data_reg;
clear data_reg

data = double(data);


% %cells 
% 
% CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
% cd(CD);
% %run FlahsingStim_dataStruct.m to organize data first
% load('dataMasks.mat')
% mask_cell = dataMasks.retTuning;
% 
% dataTC = stackGetTimeCourses(data,mask_cell);


    % variables
nTrials = input.trialSinceReset;
nTrials = nTrials-1;
cLeverDown = double(cell2mat(input.cLeverDown));
cLeverDown = cLeverDown(1:end-1);
cTargetOn = input.cTargetOn;
for itrial = 1:nTrials
    if isempty(cTargetOn{itrial})
        cTargetOn{itrial} = NaN;
    end
end
cTargetOn = (double(cell2mat_padded(cTargetOn)))'; %For now NaNs == 0, may need to change...
cTargetOn = cTargetOn(1:end-1);
cLeverUp = double(cell2mat(input.cLeverUp));
cLeverUp = cLeverUp(1:end-1);
tCyclesOn = double(cell2mat(input.tCyclesOn));
tCyclesOn = tCyclesOn(1:end-1);
% ONms = input.stimOnTimeMs;
% OFFms = input.stimOffTimeMs;
ONfr = input.nFramesOn;
OFFfr = input.nFramesOff;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
Block2ON = Block2ON(1:end-1);
TrialOutcome = input.trialOutcomeCell;
TrialOutcome = TrialOutcome(1:end-1);
Cycles = unique(tCyclesOn);
minCyclesOn = input.minCyclesOn;
maxCyclesOn = input.maxCyclesOn;
absLeverDown = 5;
absTargetOn = cTargetOn - cLeverDown;
for itrial = 1:nTrials
    if absTargetOn(itrial) < 0
        absTargetOn(itrial) = NaN;
    end
end    
absLeverUp = cLeverUp - cLeverDown;

for itrial = 1:nTrials
    trialInd(itrial,:) = cLeverDown(itrial)-30:cLeverDown(itrial)+29;
end

for itrial = 1:nTrials
    trialInd_F(itrial,:) = cLeverDown(itrial)-30:cLeverDown(itrial);
end

start = 1;
for itrial = 1:nTrials
    dataF(:,:,(start:start+59)) = data(:,:,trialInd(itrial,:));
    start = start +60;
end

for itrial = 1:nTrials
    F(:,:,itrial) = mean(data(:,:,trialInd_F(itrial,:)),3);
end

start = 1;
for itrial = 1:nTrials
    dF(:,:,(start:start+59)) = bsxfun(@minus,data(:,:,trialInd(itrial,:)),F(:,:,itrial));
    start = start+60;
end

start = 1;
for itrial = 1:nTrials
    dFoverF(:,:,start:start+59) = bsxfun(@rdivide,dF(:,:,start:start+59),F(:,:,itrial));
    start = start+60;
end

b = 5;
siz = size(data);
corr_map = zeros(siz(1),siz(2));
for ix = b:siz(2)-b
    for iy = b:siz(1)-b
        TC = data(iy,ix,:);
        surround = (data(iy-1,ix-1,:)+data(iy-1,ix,:)+data(iy-1,ix+1,:)+data(iy,ix-1,:)+data(iy,ix+1,:)+data(iy+1,ix-1,:)+data(iy+1,ix,:)+data(iy+1,ix+1,:))/8;
        R = corrcoef(TC,surround);
        corr_map(iy,ix) = R(1,2);
    end
end

figure; imagesq(corr_map); colormap(gray)
bwout = imCellEditInteractive(corr_map);
mask_cell = bwlabel(bwout);

CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
%run FlahsingStim_dataStruct.m to organize data first
load('dataMasks.mat')
mask_cell = dataMasks.dirTuning;

dataTC = stackGetTimeCourses(dFoverF,mask_cell);
dataTC_dF = stackGetTimeCourses(dF,mask_cell);
dataTC_F = stackGetTimeCourses(dataF,mask_cell);


start = 1;
for itrial = 1:nTrials
    dataCellsTrials(:,:,itrial) = dataTC(start:start+59,:);
    start = start+60;
end

% for itrial = 1:nTrials
%     trialInd(itrial,:) = cLeverDown(itrial)-30:cLeverDown(itrial)+29;
% end
% 
% for itrial = 1:nTrials
%     F(itrial,:) = mean(dataTC(trialInd(:,itrial),:),1);
% end
% 
% start = 1;
% for itrial = 1:nTrials
%     dF(start:start+59,:) = bsxfun(@minus,dataTC((trialInd(itrial,:))'),F(itrial,:));
%     start = start+60;
% end
% 
% start = 1;
% for itrial = 1:nTrials
%     dFoverF(start:start+59,:) = bsxfun(@rdivide,dF(start:start+59,:),F(itrial,:));
%     start = start+60;
% end
% 
% start = 1;
% for itrial = 1:nTrials
%     dataCellsTrials(:,:,itrial) = dFoverF(start:start+59,:);
%     start = start+60;
% end

%sort auditory and visual trials
A_ind = find(Block2ON == 1);
V_ind = find(Block2ON ==0);

dataCells_meanTrials_A = mean(dataCellsTrials(:,:,A_ind),3);
dataCells_meanTrials_V = mean(dataCellsTrials(:,:,V_ind),3);

data_meanCellsTrials_A = mean(dataCells_meanTrials_A,2);
data_meanCellsTrials_V = mean(dataCells_meanTrials_V,2);

figure;
plot(data_meanCellsTrials_A,'r');
hold on
plot(data_meanCellsTrials_V,'g');
hold on
vline(30,'k');

start = 1;
figure;
for iplot = 1:16
    subplot(4,4,iplot);
    plot(dataCells_meanTrials_A(:,start),'r');
    hold on
    plot(dataCells_meanTrials_V(:,start),'g');
    hold on
    vline(30,'k');
    start = start+1;
end


