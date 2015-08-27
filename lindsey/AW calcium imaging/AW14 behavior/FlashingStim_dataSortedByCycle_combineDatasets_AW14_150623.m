%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
SubNum = '614';
mouse = 'AW14';
date = '150623';
time_mat = ['1217'];
runs = ['004'];
ImgFolder = '004';
nrun = size(runs,1);
frame_rate = 30;

runstr = runs(1,:);
if nrun>1
    for irun = 2:nrun
        runstr = [runstr '-' runs(irun,:)];
    end
end
fnout = ['Z:\home\lindsey\Analysis\Behavior\EyeTracking\' mouse '-' date '\' mouse '-' date '-' runstr];
%% load and combine mworks data
input = [];
for irun = 1:nrun
    time = time_mat(irun,:);
    fn_mworks = ['\\CRASH.dhe.duke.edu\data\home\andrew\Behavior\Data\data-i' SubNum '-' date '-' time '.mat'];
    if irun == 1
        input = mwLoadData(fn_mworks, [], []);
    else
        input = [input mwLoadData(fn_mworks, [], [])];
    end
end
input = concatenateDataBlocks(input);

run_trials = input.trialsSinceReset;
cLeverDown = cell2mat(input.cLeverDown);
cLeverUp = cell2mat(input.cLeverUp);
cTargetOn = celleqel2mat_padded(input.cTargetOn);
cItiStart = cell2mat(input.cItiStart);

for irun = 1:nrun
    if irun < nrun
        offset = size(Area{irun},1);
        startTrial = run_trials(irun)+1;
        endTrial = run_trials(irun)+run_trials(irun+1);
        cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
        cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
        cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
        cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
    end
end
ntrials = length(input.trialOutcomeCell);
trialOutcome = cell2mat(input.trialOutcomeCell);
tCyclesOn = cell2mat(input.tCyclesOn);
cycles = unique(tCyclesOn);
V_ind = find(cell2mat(input.tBlock2TrialNumber) == 0);
AV_ind = find(cell2mat(input.tBlock2TrialNumber) == 1);
cycTime = input.nFramesOn+input.nFramesOff;

%% load and combine timecourse
dataTC = [];
for irun = 1:nrun
    fnTC = fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fnTC);
    % load('dataTC.mat');
    load('Timecourses.mat')
    dataTC = cat(2, dataTC, dataTimecourse.dataTCsub);
end

clear dataTimecourse

%% divide up data by cycle- align to lever down

for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));
    Data = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDFoverF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    if cLeverDown(end,1)+30 > size(dataTC,1)
        for itrial = 1:length(ind)-1
            Data(:,:,itrial) = dataTC(cLeverDown(ind(itrial))-30:cLeverDown(ind(itrial))+29+(cycTime.*(cycles(icyc)+1)),:);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
        end
    else
        for itrial = 1:length(ind)
            Data(:,:,itrial) = dataTC(cLeverDown(ind(itrial))-30:cLeverDown(ind(itrial))+29+(cycTime.*(cycles(icyc)+1)),:);
            DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
            DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
        end
    end
    cycData{icyc} = Data;
    cycDataDF{icyc} = DataDF;
    cycDataDFoverF{icyc} = DataDFoverF;
end

%% Align data to lever up
Data = zeros(105,size(dataTC,2),ntrials);
DataDF = zeros(105,size(dataTC,2),ntrials);
DataDFoverF = zeros(105,size(dataTC,2),ntrials);
if cLeverUp(end,1)+30 > size(dataTC,1)
    for itrial = 1:ntrials-1
        Data(:,:,itrial) = dataTC(cLeverUp(itrial)-30:cLeverUp(itrial)+74,:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
else
    for itrial = 1:ntrials
        Data(:,:,itrial) = dataTC(cLeverUp(itrial)-30:cLeverUp(itrial)+74,:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
end

DataDFoverFavg = squeeze(mean(DataDFoverF,2));
FIx = find(strcmp(input.trialOutcomeCell, 'failure'));
SIx = find(strcmp(input.trialOutcomeCell, 'success'));
FIxlong = intersect(find(tCyclesOn>3), FIx);
SIxlong = intersect(find(tCyclesOn>3), SIx);
Fb1Ix = intersect(V_ind, FIx);
Fb2Ix = intersect(AV_ind, FIx);
Sb1Ix = intersect(V_ind, SIx);
Sb2Ix = intersect(AV_ind, SIx);

figure;
tt = [-30:74].*(1000/30);
subplot(2,2,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k'); 
hold on
vline(0,'k')
title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
xlim([-200 500])
subplot(2,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
hold on
vline(0,'k')
title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
xlim([-200 500])
subplot(2,2,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'g'); 
hold on
vline(0,'k')
title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
xlim([-200 500])
subplot(2,2,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'g'); 
hold on
vline(0,'k')
alignYaxes
title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
xlim([-200 500])



   
