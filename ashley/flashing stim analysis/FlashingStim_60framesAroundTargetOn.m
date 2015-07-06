SubNum = '604';
mouse = 'AW04';
date = '141022';
time = '1600';
ImgFolder = '002';

% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

%load timecourse
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('dataTC.mat');
dataTC = dataTimecourse.dataTC;

%variables from mworks
cLeverDown = cell2mat_padded(input.cLeverDown);
cTargetOn = cell2mat_padded(input.cTargetOn);
tCyclesOn = cell2mat_padded(input.tCyclesOn);
block2 = cell2mat_padded(input.tBlock2TrialNumber);
cycles = unique(tCyclesOn);
% cycTime = input.nFramesOn + input.nFramesOff;
frameRateS = 30; %hard-coded for now, but should be available in scanbox-yeti datasets' info file
RateFRperMS = frameRateS/1000;
cycTime = ceil((input.stimOnTimeMs+input.stimOffTimeMs)*RateFRperMS);
nTrials = input.trialSinceReset;

% special case where last trial take place within last 2 frames collected
% by mworks, but not scanbox
cLeverDown = cLeverDown(1:end-10,:);
cTargetOn = cTargetOn(1:end-10,:);
tCyclesOn = tCyclesOn(1:end-10,:);
block2 = block2(1:end-10,:);
nTrials = nTrials-10;

%dFoverF
F = zeros(30,size(dataTC,2),nTrials);
Data = zeros(60,size(dataTC,2),nTrials);
DataDF = zeros(60,size(dataTC,2),nTrials);
DFoverF = zeros(60,size(dataTC,2),nTrials);
for itrial = 1:nTrials
    F(:,:,itrial) = dataTC(cLeverDown(itrial)-30:cLeverDown(itrial)-1,:);
    Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+29,:);
    DataDF(:,:,itrial) = bsxfun(@minus,Data(:,:,itrial),mean(F(:,:,itrial),1));
    DFoverF(:,:,itrial) = bsxfun(@rdivide,DataDF(:,:,itrial),mean(F(:,:,itrial),1));
end

% save('DFoverF_60framesTargetOn.mat','DFoverF');

%plot averages
visualIndex = find(block2==0);
auditoryIndex = find(block2==1);
DFoverFmeanV = mean(mean(DFoverF(:,:,visualIndex),3),2);
DFoverFmeanA = mean(mean(DFoverF(:,:,auditoryIndex),3),2);
figure;
plot(DFoverFmeanV,'g');
hold on
plot(DFoverFmeanA,'r');
hold on
vline(30,'k');

%plot average for each cycle length
DFoverF_cycAvgV = zeros(60,length(cycles));
DFoverF_cycAvgA = zeros(60,length(cycles));
for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));
    v_ind = intersect((find(block2 == 0)),ind);
    a_ind = intersect((find(block2 == 1)),ind);
    DFoverF_cycAvgV(:,icyc) = mean(mean(DFoverF(:,:,v_ind),3),2);
    DFoverF_cycAvgA(:,icyc) = mean(mean(DFoverF(:,:,a_ind),3),2);
end


CV = colormap(brewermap(2*length(cycles),'Greens'));
CA = colormap(brewermap(2*length(cycles),'YlOrRd'));
figure;
for icyc = 1:length(cycles)
    plot(DFoverF_cycAvgV(:,icyc),'color',CV(icyc+length(cycles),:))
    hold on
    plot(DFoverF_cycAvgA(:,icyc),'color',CA(icyc+length(cycles),:))
    hold on
    vline(30,'k')
    hold on
end

figure;
plot(DFoverF_cycAvgV(:,1),'color',CV(1+length(cycles),:))
hold on
plot(DFoverF_cycAvgV(:,length(cycles)),'color',CV(length(cycles)+length(cycles),:));
hold on
plot(DFoverF_cycAvgA(:,1),'color',CA(1+length(cycles),:))
hold on
plot(DFoverF_cycAvgA(:,length(cycles)),'color',CA(length(cycles)+length(cycles),:));
hold on
vline(30,'k')

%save
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
save('DFoverF_60framesTargetOn.mat','DFoverF');