%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
SubNum = '613';
mouse = 'AW13';
date = '150511';
time_mat = ['1420'; '1437'; '1453'];
runs = ['001'; '002'; '003'];
nrun = size(runs,1);
frame_rate = 30;

runstr = runs(1,:);
if nrun>1
    for irun = 2:nrun
        runstr = [runstr '-' runs(irun,:)];
    end
end
%% load and combine mworks data and timecourses
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

dataTC = [];
for irun = 1:nrun
    ImgFolder = runs(irun,:);
    fnTC = fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fnTC);
    load('Timecourses.mat')
    dataTC = cat(2, dataTC, dataTimecourse.dataTCsub);
    if irun < nrun
        offset = size(dataTimecourse.dataTCsub,3);
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

tGratingDirectionDeg = chop(cell2mat(input.tGratingDirectionDeg),4);
Dirs = unique(tGratingDirectionDeg);
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
MIx = find(strcmp(input.trialOutcomeCell, 'ignore'));
FIxlong = intersect(find(tCyclesOn>3), FIx);
SIxlong = intersect(find(tCyclesOn>3), SIx);
MIxlong = intersect(find(tCyclesOn>3), MIx);
Fb1Ix = intersect(V_ind, FIx);
Fb2Ix = intersect(AV_ind, FIx);
Sb1Ix = intersect(V_ind, SIx);
Sb2Ix = intersect(AV_ind, SIx);
Mb1Ix = intersect(V_ind, MIx);
Mb2Ix = intersect(AV_ind, MIx);

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


%% Align data to target on up
Data = zeros(105,size(dataTC,2),ntrials);
DataDF = zeros(105,size(dataTC,2),ntrials);
DataDFoverF = zeros(105,size(dataTC,2),ntrials);

for itrial = 1:ntrials-1
    if ~isnan(cTargetOn(itrial))
        Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+74,:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
end

DataDFoverFavg = squeeze(mean(DataDFoverF,2));
figure;
for idir = 1:length(Dirs)
    subplot(3,3,idir)
    ind1 = intersect(Sb1Ix, find(tGratingDirectionDeg==Dirs(idir)));
    ind2 = intersect(Mb1Ix, find(tGratingDirectionDeg==Dirs(idir)));
    plot(tt, nanmean(DataDFoverFavg(:,ind1),2), '-k');
    hold on;
    plot(tt, nanmean(DataDFoverFavg(:,ind2),2), '-r');
    ylim([-0.02 0.04])
    hold on
    vline(0,'k')
    hold on
    vline(550,'k')
    title([num2str(length(ind1)) ' - ' num2str(length(ind2))])
end

   
