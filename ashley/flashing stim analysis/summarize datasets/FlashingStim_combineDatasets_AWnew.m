% choose dataset to run
% close all
% clear all
awFSAVdatasets
try
    exp
catch
    error('choose dataset from awFSAVdatasets and set exp variable')
end
SubNum = expt(exp).SubNum;
mouse = expt(exp).mouse;
date = expt(exp).date;
runs = expt(exp).runs;
nrun = size(runs,1);
expFolder = expt(exp).folder;
time_mat = expt(exp).time_mat;

%%
%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
runstr = runs(1,:);
if nrun>1
    for irun = 2:nrun
        runstr = [runstr '-' runs(irun,:)];
    end
end
% fnout = fullfile('Z:\home\lindsey\Analysis\2P', mouse, date, [date '_' mouse '_' runstr '_']);
fnout = ['Z:\Analysis\' mouse '\' expFolder '\' date '\' mouse '-' date '-' runstr];

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

%%
run_trials = input.trialsSinceReset;
cLeverDown = cell2mat(input.cLeverDown);
cLeverUp = cell2mat(input.cLeverUp);
cTargetOn = celleqel2mat_padded(input.cTargetOn);
cStimOn = celleqel2mat_padded(input.cStimOn);
cItiStart = cell2mat(input.cItiStart);
if expt(exp).catch == 1;
    cCatchOn = celleqel2mat_padded(input.cCatchOn);
    isFA = celleqel2mat_padded(input.tFalseAlarm);
end

dataTC = [];
offset = 0;
for irun = 1:nrun
    ImgFolder = runs(irun,:);
    fnTC = fullfile('\\CRASH.dhe.duke.edu\data\home\ashley\analysis\',mouse,'two-photon imaging', date, ImgFolder);
    cd(fnTC);
    load('Timecourses.mat')
    dataTC = cat(1, dataTC, dataTimecourse.dataTCsub);
    offset = offset+size(dataTimecourse.dataTCsub,1);
    if irun < nrun
        startTrial = sum(run_trials(1, 1:irun),2)+1;
        endTrial = sum(run_trials(1,1:irun+1),2);
        cLeverDown(1,startTrial:endTrial) = cLeverDown(1,startTrial:endTrial)+offset;
        cLeverUp(1,startTrial:endTrial) = cLeverUp(1,startTrial:endTrial)+offset;
        cTargetOn(1,startTrial:endTrial) = cTargetOn(1,startTrial:endTrial)+offset;
        cStimOn(1,startTrial:endTrial) = cStimOn(1,startTrial:endTrial)+offset;
        cItiStart(1,startTrial:endTrial) = cItiStart(1,startTrial:endTrial)+offset;
        if expt(exp).catch == 1;
            cCatchOn(1,startTrial:endTrial) = cCatchOn(1,startTrial:endTrial)+offset;
            isFA(1,startTrial:endTrial) = cCatchOn(1,startTrial:endTrial)+offset;
        end
    end
end

ntrials = length(input.trialOutcomeCell);
trialOutcome = cell2mat(input.trialOutcomeCell);
tCyclesOn = cell2mat(input.tCyclesOn);
nCyclesOn = cell2mat(input.nCyclesOn);

cycles = unique(tCyclesOn);
V_ind = find(cell2mat(input.tBlock2TrialNumber) == 0);
AV_ind = find(cell2mat(input.tBlock2TrialNumber) == 1);
cycTime = input.nFramesOn+input.nFramesOff;
tooFastTime = input.nFramesTooFast;
maxReactTime = input.nFramesReact;

DirectionDeg = cell2mat(input.tGratingDirectionDeg);
Dirs = unique(DirectionDeg);
tGratingDirectionDeg = chop(cell2mat(input.tGratingDirectionDeg),4);

if expt(exp).catch == 1;
    catchCycle = cell2mat(input.catchCyclesOn);
    catchDirectionDeg = cell2mat_padded(input.tCatchGratingDirectionDeg);
    catchDirs = unique(catchDirectionDeg);
    isCatchTrial = catchDirectionDeg > 0;
    catchTrialOutcome = num2cell(NaN(length(nCyclesOn),1));
    catchIndex = find(isCatchTrial == 1);
    for i = 1:sum(isCatchTrial)
        if isFA(catchIndex(i)) == 1
            catchTrialOutcome{catchIndex(i),1} = 'FA';
        elseif cCatchOn(catchIndex(i)) == 0
            catchTrialOutcome{catchIndex(i),1} = 'failure';
        elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) < tooFastTime
            catchTrialOutcome{catchIndex(i),1} = 'failure';
        elseif (cLeverUp(catchIndex(i)) - cCatchOn(catchIndex(i))) > maxReactTime
            catchTrialOutcome{catchIndex(i),1} = 'CR';
        end
    end
end

clear dataTimecourse

try
    cd([fnout '\' runstr])
catch
    mkdir(fnout,runstr)
end

clear date
save([mouse '-' expt(exp).date '-' runstr '-' comboInputDataTCplusVar]);
save([date mouse '-' expt(exp).date '-' runstr '-' comboInputDataTCplusVar]);

%% divide up data by cycle- align to lever down

for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));
    Data = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    DataDFoverF = zeros(cycTime.*(cycles(icyc)+1)+60,size(dataTC,2),length(ind));
    if cLeverDown(end)+30 > size(dataTC,1)
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

DataDFoverFavg = squeeze(mean(DataDFoverF(:,:,:),2));
FIx = find(strcmp(input.trialOutcomeCell, 'failure'));
SIx = find(strcmp(input.trialOutcomeCell, 'success'));
MIx = find(strcmp(input.trialOutcomeCell, 'ignore'));
FIxlong = intersect(find(tCyclesOn>3), FIx);
SIxlong = intersect(find(tCyclesOn>3), SIx);
MIxlong = intersect(find(tCyclesOn>3), MIx);
Fb1Ix = intersect(V_ind, FIxlong);
Fb2Ix = intersect(AV_ind, FIxlong);
Sb1Ix = intersect(V_ind, SIxlong);
Sb2Ix = intersect(AV_ind, SIxlong);
Mb1Ix = intersect(V_ind, MIxlong);
Mb2Ix = intersect(AV_ind, MIxlong);
Rb1Ix = intersect(V_ind, find(tCyclesOn>3));
Rb2Ix = intersect(AV_ind, find(tCyclesOn>3));

figure;
tt = [-30:74].*(1000/30);
subplot(2,2,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
subplot(2,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
subplot(2,2,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
subplot(2,2,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
alignYaxes
title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
suptitle([date ' ' mouse ' ' runstr '- align to lever release- selective cells'])
print([fnout 'release_align_SE_AV_pref.pdf'], '-dpdf')

%% Align data to previous baseline stim on
Data = zeros(105,size(dataTC,2),ntrials);
DataDF = zeros(105,size(dataTC,2),ntrials);
DataDFoverF = zeros(105,size(dataTC,2),ntrials);
if cStimOn(end,1)+30 > size(dataTC,1)
    for itrial = 1:ntrials-1
        Data(:,:,itrial) = dataTC(cStimOn(itrial)-30:cStimOn(itrial)+74,:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
else
    for itrial = 1:ntrials
        Data(:,:,itrial) = dataTC(cStimOn(itrial)-30:cStimOn(itrial)+74,:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
end

DataDFoverFavg = squeeze(mean(DataDFoverF,2));
figure;
tt = [-30:74].*(1000/30);
subplot(2,2,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
subplot(2,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
subplot(2,2,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
subplot(2,2,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
alignYaxes
title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
suptitle([date ' ' mouse ' ' runstr '- align to previous baseline stim on'])
print([fnout 'baselineStim_align_SE_AV.pdf'], '-dpdf')

%% Align data to previous stim on
Data = zeros(105,size(dataTC,2),ntrials);
Data_CR = zeros(105,size(dataTC,2),ntrials);
DataDF = zeros(105,size(dataTC,2),ntrials);
DataDFoverF = zeros(105,size(dataTC,2),ntrials);
DataDF_CR = zeros(105,size(dataTC,2),ntrials);
DataDFoverF_CR = zeros(105,size(dataTC,2),ntrials);
if cTargetOn(end,1)+30 > size(dataTC,1)
    for itrial = 1:ntrials-1
        if strcmp(input.trialOutcomeCell(itrial), 'failure')
            Data(:,:,itrial) = dataTC(cStimOn(itrial)-30:cStimOn(itrial)+74,:);
            Data_CR(:,:,itrial) = dataTC(cStimOn(itrial)-30-11:cStimOn(itrial)+74-11,:);
        else
            Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+74,:);
            Data_CR(:,:,itrial) = dataTC(cTargetOn(itrial)-30-11:cTargetOn(itrial)+74-11,:);
        end
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDF_CR(:,:,itrial) = bsxfun(@minus, Data_CR(:,:,itrial), mean(Data_CR(1:30,:,itrial),1));
        DataDFoverF_CR(:,:,itrial) = bsxfun(@rdivide, DataDF_CR(:,:,itrial), mean(Data_CR(1:30,:,itrial),1));
    end
else
    for itrial = 1:ntrials
        if strcmp(input.trialOutcomeCell(itrial), 'failure')
            Data(:,:,itrial) = dataTC(cStimOn(itrial)-30:cStimOn(itrial)+74,:);
            Data_CR(:,:,itrial) = dataTC(cStimOn(itrial)-30-11:cStimOn(itrial)+74-11,:);
        else
            Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+74,:);
            Data_CR(:,:,itrial) = dataTC(cTargetOn(itrial)-30-11:cTargetOn(itrial)+74-11,:);
        end
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDF_CR(:,:,itrial) = bsxfun(@minus, Data_CR(:,:,itrial), mean(Data_CR(1:30,:,itrial),1));
        DataDFoverF_CR(:,:,itrial) = bsxfun(@rdivide, DataDF_CR(:,:,itrial), mean(Data_CR(1:30,:,itrial),1));
    end
end

for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
figure;
tt = [-30:74].*(1000/30);
subplot(2,2,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
subplot(2,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
subplot(2,2,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
subplot(2,2,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
alignYaxes
title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
suptitle([date ' ' mouse ' ' runstr '- align to previous stim on- ' num2str(Oris(iOri)) ' deg select: n = ' num2str(length(cellsSelect{iOri}))])
print([fnout 'prevStimOn_align_SE_AV_' num2str(Oris(iOri)) 'pref.pdf'], '-dpdf')
end

colmat = strvcat('k', 'g', 'r', 'b');
OriNs = [];
figure;
for iOri = 1:nOri
OriNs = [OriNs length(cellsSelect{iOri})];
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,1)
A(iOri) = shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), colmat(iOri,:));
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual Hit trials: n = ' num2str(length(Sb1Ix))]) 
subplot(2,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)),  colmat(iOri,:));
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual FA trials: n = ' num2str(length(Fb1Ix))]) 
subplot(2,2,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), colmat(iOri,:));
xlim([-500 1500])
hold on
vline(0,'k')
title(['Auditory Hit trials: n = ' num2str(length(Sb2Ix))]) 
subplot(2,2,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), colmat(iOri,:));
xlim([-500 1500])
hold on
vline(0,'k')
alignYaxes
title(['Auditory FA trials: n = ' num2str(length(Fb2Ix))]) 
end
legend([A(1).mainLine, A(2).mainLine, A(3).mainLine, A(4).mainLine], num2str(Oris'))
suptitle(sprintf([date ' ' mouse ' ' runstr '- align to previous stim on- \n' num2str(Oris) ' deg select: n = ' num2str(OriNs)]))
print([fnout 'prevStimOn_align_SE_AV_eachOripref.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,1)
plot(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), colmat(iOri,:));
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual Hit trials: n = ' num2str(length(Sb1Ix))]) 
subplot(2,2,2)
plot(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2),  colmat(iOri,:));
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual FA trials: n = ' num2str(length(Fb1Ix))]) 
subplot(2,2,3)
plot(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), colmat(iOri,:));
xlim([-500 1500])
hold on
vline(0,'k')
title(['Auditory Hit trials: n = ' num2str(length(Sb2Ix))]) 
subplot(2,2,4)
plot(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), colmat(iOri,:));
xlim([-500 1500])
hold on
vline(0,'k')
title(['Auditory FA trials: n = ' num2str(length(Fb2Ix))]) 
end
alignYaxes
legend(num2str(Oris'))
suptitle(sprintf([date ' ' mouse ' ' runstr '- align to previous stim on- \n' num2str(Oris) ' deg select: n = ' num2str(OriNs)]))
print([fnout 'prevStimOn_align_SE_AV_eachOripref_noerror.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
DataDFoverFavg_CR = squeeze(mean(DataDFoverF_CR(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,1)
plot(tt, nanmean(DataDFoverFavg(:,Mb1Ix),2), colmat(iOri,:));
xlim([-500 1500])
ylim([-0.01 0.05])
hold on
vline(0,'k')
title(['Visual Miss trials: n = ' num2str(length(Mb1Ix))]) 
subplot(2,2,2)
plot(tt, nanmean(DataDFoverFavg_CR(:,Rb1Ix),2),  colmat(iOri,:));
xlim([-500 1500])
ylim([-0.01 0.05])
hold on
vline(0,'k')
title(['Visual CR trials: n = ' num2str(length(Rb1Ix))]) 
subplot(2,2,3)
plot(tt, nanmean(DataDFoverFavg(:,Mb2Ix),2), colmat(iOri,:));
xlim([-500 1500])
ylim([-0.01 0.05])
hold on
vline(0,'k')
title(['Auditory Miss trials: n = ' num2str(length(Mb2Ix))]) 
subplot(2,2,4)
plot(tt, nanmean(DataDFoverFavg_CR(:,Rb2Ix),2), colmat(iOri,:));
xlim([-500 1500])
ylim([-0.01 0.05])
hold on
vline(0,'k')
title(['Auditory CR trials: n = ' num2str(length(Rb2Ix))]) 
end
alignYaxes
legend(num2str(Oris'))
suptitle(sprintf([date ' ' mouse ' ' runstr '- align to previous stim on- \n' num2str(Oris) ' deg select: n = ' num2str(OriNs)]))
print([fnout 'prevStimOn_align_CM_AV_eachOripref_noerror.pdf'], '-dpdf')


targs = unique(tGratingDirectionDeg);
nTarg = length(targs);
n = ceil(sqrt(nTarg-1));
if (n^2)-n<nTarg-1
    n2 = n;
else
    n2 = n-1;
end
figure;
for iOri = 1:nOri
    DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
    tt = [-30:74].*(1000/30);
    for itarg = 2:nTarg
        Sb1Ixtarg = intersect(Sb1Ix, find(tGratingDirectionDeg == targs(itarg)));
        subplot(n,n2,itarg-1)
        plot(tt, nanmean(DataDFoverFavg(:,Sb1Ixtarg),2), colmat(iOri,:));
        xlim([-500 1500])
        ylim([-.02 .05])
        hold on
        vline(0,'k')
        title([num2str(targs(itarg)) ' deg trials: n = ' num2str(length(Sb1Ixtarg))]) 
    end
end
legend(num2str(Oris'))
suptitle(sprintf([date ' ' mouse ' ' runstr '- Visual Hits: align to previous stim on- \n' num2str(Oris) ' deg select: n = ' num2str(OriNs)]))
print([fnout 'prevStimOn_align_Hits_eachOripref_noerror.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
    DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
    tt = [-30:74].*(1000/30);
    for itarg = 2:nTarg
        Mb1Ixtarg = intersect(Mb1Ix, find(tGratingDirectionDeg == targs(itarg)));
        subplot(n,n2,itarg-1)
        plot(tt, nanmean(DataDFoverFavg(:,Mb1Ixtarg),2), colmat(iOri,:));
        xlim([-500 1500])
        ylim([-.02 .05])
        hold on
        vline(0,'k')
        title([num2str(targs(itarg)) ' deg trials: n = ' num2str(length(Mb1Ixtarg))]) 
    end
end
legend(num2str(Oris'))
suptitle(sprintf([date ' ' mouse ' ' runstr '- Visual Misses: align to previous stim on- \n' num2str(Oris) ' deg select: n = ' num2str(OriNs)]))
print([fnout 'prevStimOn_align_Misses_eachOripref_noerror.pdf'], '-dpdf')

tt = [-30:74].*(1000/30);
for iOri = 1:nOri
    figure;
    for itarg = 2:nTarg
        Sb1Ixtarg = intersect(Sb1Ix, find(tGratingDirectionDeg == targs(itarg)));
        Mb1Ixtarg = intersect(Mb1Ix, find(tGratingDirectionDeg == targs(itarg)));
        subplot(n,n2,itarg-1)
        DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
        if length(Sb1Ixtarg)>2
            shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ixtarg),2), nanstd(DataDFoverFavg(:,Sb1Ixtarg),[],2)/sqrt(length(Sb1Ixtarg)), 'k');
        else
            plot(tt, nanmean(DataDFoverFavg(:,Sb1Ixtarg),2), 'k');
        end
        hold on
        if length(Mb1Ixtarg)>2
            shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Mb1Ixtarg),2), nanstd(DataDFoverFavg(:,Mb1Ixtarg),[],2)/sqrt(length(Mb1Ixtarg)), 'r');
        else
            plot(tt, nanmean(DataDFoverFavg(:,Mb1Ixtarg),2), 'r');
        end
        xlim([-500 1500])
        ylim([-.02 .05])
        title([num2str(targs(itarg)) ' deg trials: Hit n = ' num2str(length(Sb1Ixtarg)) '; Miss n = ' num2str(length(Mb1Ixtarg)) ]) 
    end
    suptitle(sprintf([date ' ' mouse ' ' runstr '- Visual Trials: align to previous stim on- \n ' num2str(Oris(iOri)) ' deg Selective: n = ' num2str(length(cellsSelect{iOri}))]))
    print([fnout 'prevStimOn_align_HM_ByTarget_' num2str(Oris(iOri)) 'deg.pdf'], '-dpdf')
end

tt = [-30:74].*(1000/30);
figure;
kmat = flipud(repmat([0.2:0.1:1],[3, 1])');
for iOri = 1:nOri
    subplot(2,2,iOri)
    DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
    shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
    hold on
    Sb1IxtargN = [];
    for itarg = 2:nTarg
        Sb1Ixtarg = intersect(Sb1Ix, find(tGratingDirectionDeg == targs(itarg)));
        Sb1IxtargN = [Sb1IxtargN length(Sb1Ixtarg)];
        B(itarg) = plot(tt, nanmean(DataDFoverFavg(:,Sb1Ixtarg),2), 'Color', kmat(itarg,:));
        hold on
        xlim([-500 1500])
        ylim([-.02 .05])
    end
    title([num2str(Oris(iOri)) ' Pref cells- n = ' num2str(OriNs(iOri))]) 
end
legend(B(2:nTarg), [num2str(chop(targs(2:nTarg),2)') repmat(' n =', [nTarg-1, 1]) num2str(Sb1IxtargN')])
suptitle(sprintf([date ' ' mouse ' ' runstr '- Visual Trials: align to previous stim on- \n Hits: ' num2str(length(Sb1Ix)) '; FAs n = ' num2str(length(Fb1Ix))]))
print([fnout 'prevStimOn_align_ByOri_ByTarget_HF.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
    subplot(2,2,iOri)
    DataDFoverFavg_CR = squeeze(mean(DataDFoverF_CR(:,cellsSelect{iOri},:),2));
    shadedErrorBar(tt, nanmean(DataDFoverFavg_CR(:,Rb1Ix),2), nanstd(DataDFoverFavg_CR(:,Rb1Ix),[],2)/sqrt(length(Rb1Ix)), 'r');
    hold on
    Mb1IxtargN = [];
    DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
    for itarg = 2:nTarg
        Mb1Ixtarg = intersect(Mb1Ix, find(tGratingDirectionDeg == targs(itarg)));
        Mb1IxtargN = [Mb1IxtargN length(Mb1Ixtarg)];
        B(itarg) = plot(tt, nanmean(DataDFoverFavg(:,Mb1Ixtarg),2), 'Color', kmat(itarg,:));
        hold on
        xlim([-500 1500])
        ylim([-.02 .05])
    end
    title([num2str(Oris(iOri)) ' Pref cells- n = ' num2str(OriNs(iOri))]) 
end
legend(B(2:nTarg), [num2str(chop(targs(2:nTarg),2)') repmat(' n =', [nTarg-1, 1]) num2str(Mb1IxtargN')])
suptitle(sprintf([date ' ' mouse ' ' runstr '- Visual Trials: align to previous stim on- \n Miss: ' num2str(length(Mb1Ix)) '; CRs n = ' num2str(length(Rb1Ix))]))
print([fnout 'prevStimOn_align_ByOri_ByTarget_MC.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(1,2,1)
plot(tt, nanmean(DataDFoverFavg(:,Mb1Ix),2), colmat(iOri,:));
xlim([-500 1500])
ylim([-0.01 0.03])
hold on
vline(0,'k')
title(['Visual Miss trials: n = ' num2str(length(Mb1Ix))])
subplot(1,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Mb1Ix),2), nanstd(DataDFoverFavg(:,Mb1Ix),[],2)/sqrt(length(Mb1Ix)), colmat(iOri,:));
xlim([-500 1500])
ylim([-0.01 0.03])
hold on
vline(0,'k')
title(['Visual Miss trials: n = ' num2str(length(Mb1Ix))]) 
end
subplot(1,2,1)
legend(num2str(Oris'))
suptitle(sprintf([date ' ' mouse ' ' runstr '- align to previous stim on- \n' num2str(Oris) ' deg select: n = ' num2str(OriNs)]))
print([fnout 'prevStimOn_align_M_Visonly_eachOripref_noerror.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,iOri)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Mb1Ix),2), nanstd(DataDFoverFavg(:,Mb1Ix),[],2)/sqrt(length(Mb1Ix)), 'r');
xlim([-500 1500])
ylim([-0.01 0.03])
hold on
vline(0,'k')
title([num2str(Oris(iOri)) ' Pref cells- n = ' num2str(OriNs(iOri))]) 
end
suptitle(['Visual trials: Hits n = ' num2str(length(Sb1Ix)) '; Misses n = ' num2str(length(Mb1Ix))])
print([fnout 'prevStimOn_align_SM_Visonly_eachOripref.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,iOri)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
xlim([-500 1500])
ylim([-0.01 0.03])
hold on
vline(0,'k')
title([num2str(Oris(iOri)) ' Pref cells- n = ' num2str(OriNs(iOri))]) 
end
suptitle(['Visual trials: Hits n = ' num2str(length(Sb1Ix)) '; FA n = ' num2str(length(Fb1Ix))])
print([fnout 'prevStimOn_align_SF_Visonly_eachOripref.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
DataDFoverFavg_CR = squeeze(mean(DataDFoverF_CR(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,iOri)
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg_CR(:,Rb1Ix),2), nanstd(DataDFoverFavg_CR(:,Rb1Ix),[],2)/sqrt(length(Rb1Ix)), 'k');
xlim([-500 1500])
ylim([-0.01 0.03])
hold on
vline(0,'k')
title([num2str(Oris(iOri)) ' Pref cells- n = ' num2str(OriNs(iOri))]) 
end
suptitle(['Visual trials: CRs n = ' num2str(length(Rb1Ix)) '; FA n = ' num2str(length(Fb1Ix))])
print([fnout 'prevStimOn_align_CF_Visonly_eachOripref.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,iOri)
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Mb2Ix),2), nanstd(DataDFoverFavg(:,Mb2Ix),[],2)/sqrt(length(Mb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k');
xlim([-500 1500])
ylim([-0.01 0.05])
hold on
vline(0,'k')
title([num2str(Oris(iOri)) ' Pref cells- n = ' num2str(OriNs(iOri))]) 
end
suptitle(['Auditory trials: Hits n = ' num2str(length(Sb2Ix)) '; Misses n = ' num2str(length(Mb2Ix))])
print([fnout 'prevStimOn_align_SM_Audonly_eachOripref.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,iOri)
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k');
xlim([-500 1500])
ylim([-0.01 0.03])
hold on
vline(0,'k')
title([num2str(Oris(iOri)) ' Pref cells- n = ' num2str(OriNs(iOri))]) 
end
suptitle(['Auditory trials: Hits n = ' num2str(length(Sb2Ix)) '; FA n = ' num2str(length(Fb2Ix))])
print([fnout 'prevStimOn_align_SF_Audonly_eachOripref.pdf'], '-dpdf')

figure;
for iOri = 1:nOri
DataDFoverFavg = squeeze(mean(DataDFoverF(:,cellsSelect{iOri},:),2));
DataDFoverFavg_CR = squeeze(mean(DataDFoverF_CR(:,cellsSelect{iOri},:),2));
tt = [-30:74].*(1000/30);
subplot(2,2,iOri)
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg_CR(:,Rb2Ix),2), nanstd(DataDFoverFavg_CR(:,Rb2Ix),[],2)/sqrt(length(Rb2Ix)), 'k');
xlim([-500 1500])
ylim([-0.01 0.03])
hold on
vline(0,'k')
title([num2str(Oris(iOri)) ' Pref cells- n = ' num2str(OriNs(iOri))]) 
end
suptitle(['Auditory trials: CRs n = ' num2str(length(Rb2Ix)) '; FA n = ' num2str(length(Fb2Ix))])
print([fnout 'prevStimOn_align_CF_Audonly_eachOripref.pdf'], '-dpdf')


DataDFoverFavg = squeeze(mean(DataDFoverF(:,notSlctvCells,:),2));
figure;
tt = [-30:74].*(1000/30);
subplot(2,2,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Mb1Ix),2), nanstd(DataDFoverFavg(:,Mb1Ix),[],2)/sqrt(length(Mb1Ix)), 'm');
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual trials: ' num2str(length(Sb1Ix)) ' Hit; ' num2str(length(Fb1Ix)) ' FA;' num2str(length(Mb1Ix)) ' Miss']) 
subplot(2,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Hit; ' num2str(length(Fb2Ix)) ' FA']) 
subplot(2,2,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['FA trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
subplot(2,2,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
alignYaxes
title(['Hit trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
suptitle([date ' ' mouse ' ' runstr '- align to previous stim on- not selective cells n = ' num2str(length(notSlctvCells))])
print([fnout 'baselineStim_align_SE_AV_notSelect.pdf'], '-dpdf')

DataDFoverFavg = squeeze(mean(DataDFoverF(:,oriSlctvCellsAll,:),2));
figure;
tt = [-30:74].*(1000/30);
subplot(2,2,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Mb1Ix),2), nanstd(DataDFoverFavg(:,Mb1Ix),[],2)/sqrt(length(Mb1Ix)), 'm');
xlim([-500 1500])
hold on
vline(0,'k')
title(['Visual trials: ' num2str(length(Sb1Ix)) ' Hit; ' num2str(length(Fb1Ix)) ' FA;' num2str(length(Mb1Ix)) ' Miss']) 
subplot(2,2,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Hit; ' num2str(length(Fb2Ix)) ' FA']) 
subplot(2,2,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
title(['FA trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
subplot(2,2,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
xlim([-500 1500])
hold on
vline(0,'k')
alignYaxes
title(['Hit trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
suptitle([date ' ' mouse ' ' runstr '- align to previous stim on- selective cells n = ' num2str(length(oriSlctvCellsAll))])
print([fnout 'baselineStim_align_SE_AV_Select.pdf'], '-dpdf')

%% align to trial start, by cycle
start = 6;
n = ceil(sqrt(length(cycles)-start+1));
if ((n^2)-n)<(length(cycles)-start+1)
    n2 = n;
else
    n2 = n-1;
end

Data = zeros((cycTime.*(cycles(end)))+30,size(dataTC,2),ntrials);
DataDF = zeros((cycTime.*(cycles(end)))+30,size(dataTC,2),ntrials);
DataDFoverF = zeros((cycTime.*(cycles(end)))+30,size(dataTC,2),ntrials);
if cLeverDown(end,1)+(cycTime.*(cycles(end))) > size(dataTC,1)
    for itrial = 1:ntrials-1
        Data(:,:,itrial) = dataTC(cLeverDown(itrial)-29:cLeverDown(itrial)+(cycTime.*(cycles(end))),:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
else
    for itrial = 1:ntrials
        Data(:,:,itrial) = dataTC(cLeverDown(itrial)-29:cLeverDown(itrial)+(cycTime.*(cycles(end))),:);
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
end

start = 6;
figure;
for icyc = start:length(cycles)
    ind = intersect(find(tCyclesOn >= cycles(icyc)), V_ind);
    subplot(n,n2,icyc-5)
    tt = [-29:(cycTime.*(cycles(icyc)))].*(1000/30);
    for iOri = 1:nOri
        DataDFoverFavg = squeeze(mean(mean(DataDFoverF(1:30+(cycTime.*(cycles(icyc))),cellsSelect{iOri},ind),2),3));
        plot(tt, DataDFoverFavg, colmat(iOri,:));
        ylim([-0.01 0.04])
        hold on
    end
    title([num2str(icyc) ' cycles- n = ' num2str(length(ind)) ' trials'])
    start = start+1;
end
suptitle(sprintf([date ' ' mouse ' ' runstr '- Visual trials align to start- \n' num2str(Oris) ' deg select: n = ' num2str(OriNs)]))
print([fnout 'trialstart_align_Visonly_byCyc_byOri.pdf'], '-dpdf')

figure;
start = 6;
for icyc = start:length(cycles)
    ind = intersect(find(tCyclesOn >= cycles(icyc)), AV_ind);
    subplot(n,n2,icyc-5)
    tt = [-29:(cycTime.*(cycles(icyc)))].*(1000/30);
    for iOri = 1:nOri
        DataDFoverFavg = squeeze(mean(mean(DataDFoverF(1:30+(cycTime.*(cycles(icyc))),cellsSelect{iOri},ind),2),3));
        plot(tt, DataDFoverFavg, colmat(iOri,:));
        ylim([-0.01 0.04])
        hold on
    end
    title([num2str(icyc) ' cycles- n = ' num2str(length(ind)) ' trials'])
    start = start+1;
end
suptitle(sprintf([date ' ' mouse ' ' runstr '- Auditory trials align to start- \n' num2str(Oris) ' deg select: n = ' num2str(OriNs)]))
print([fnout 'trialstart_align_Audonly_byCyc_byOri.pdf'], '-dpdf')

for iOri = 1:nOri
    start = 6;
    figure;
    for icyc = start:length(cycles)
        indV = intersect(find(tCyclesOn >= cycles(icyc)), V_ind);
        indAV = intersect(find(tCyclesOn >= cycles(icyc)), AV_ind);
        subplot(n,n2,icyc-5)
        tt = [-29:(cycTime.*(cycles(icyc)))].*(1000/30);
        DataDFoverFavg = squeeze(mean(mean(DataDFoverF(1:30+(cycTime.*(cycles(icyc))),cellsSelect{iOri},indAV),2),3));
        DataDFoverFsem = squeeze(std(mean(DataDFoverF(1:30+(cycTime.*(cycles(icyc))),cellsSelect{iOri},indAV),2),[],3))./sqrt(length(ind));
        shadedErrorBar(tt, DataDFoverFavg, DataDFoverFsem, 'k');
        hold on
        DataDFoverFavg = squeeze(mean(mean(DataDFoverF(1:30+(cycTime.*(cycles(icyc))),cellsSelect{iOri},indV),2),3));
        DataDFoverFsem = squeeze(std(mean(DataDFoverF(1:30+(cycTime.*(cycles(icyc))),cellsSelect{iOri},indV),2),[],3))./sqrt(length(ind));
        shadedErrorBar(tt, DataDFoverFavg, DataDFoverFsem, 'g');
        ylim([-0.01 0.05])
        title([num2str(icyc) ' cycles- n = ' num2str(length(indV)) ' Vis trials; ' num2str(length(indAV)) ' Aud trials'])
    end
    suptitle(sprintf([date ' ' mouse ' ' runstr '- align to trial start- \n' num2str(Oris(iOri)) ' deg selective cells: n = ' num2str(OriNs(iOri))]))
    print([fnout 'trialstart_align_byCyc_' num2str(Oris(iOri)) 'deg.pdf'], '-dpdf')
end

%% Align data to previous stim on, group short, med, long
Data = zeros(105,size(dataTC,2),ntrials);
DataDF = zeros(105,size(dataTC,2),ntrials);
DataDFoverF = zeros(105,size(dataTC,2),ntrials);
if cTargetOn(end,1)+30 > size(dataTC,1)
    for itrial = 1:ntrials-1
        if strcmp(input.trialOutcomeCell(itrial), 'failure')
            Data(:,:,itrial) = dataTC(cStimOn(itrial)-30:cStimOn(itrial)+74,:);
        else
            Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+74,:);
        end
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
else
    for itrial = 1:ntrials
        if strcmp(input.trialOutcomeCell(itrial), 'failure')
            Data(:,:,itrial) = dataTC(cStimOn(itrial)-30:cStimOn(itrial)+74,:);
        else
            Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+74,:);
        end
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
end

DataDFoverFavg = squeeze(mean(DataDFoverF,2));
FIx = find(strcmp(input.trialOutcomeCell, 'failure'));
SIx = find(strcmp(input.trialOutcomeCell, 'success'));
FIxlong = intersect(find(tCyclesOn>6), FIx);
SIxlong = intersect(find(tCyclesOn>6), SIx);
FIxmed = intersect(intersect(find(tCyclesOn>3), find(tCyclesOn<7)), FIx);
SIxmed = intersect(intersect(find(tCyclesOn>3), find(tCyclesOn<7)), SIx);
FIxshort = intersect(find(tCyclesOn<4), FIx);
SIxshort = intersect(find(tCyclesOn<4), SIx);

figure;
tt = [-30:74].*(1000/30);
subplot(2,3,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxshort, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxshort, V_ind)),[],2)/sqrt(length(intersect(SIxshort, V_ind))), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxshort, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxshort, V_ind)),[],2)/sqrt(length(intersect(FIxshort, V_ind))), 'r');
xlim([-500 1500])
hold on
vline(0,'k')
title(['Short visual: ' num2str(length(intersect(SIxshort, V_ind))) ' Hit; ' num2str(length(intersect(FIxshort, V_ind))) ' FA'])

subplot(2,3,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxmed, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxmed, V_ind)),[],2)/sqrt(length(intersect(SIxmed, V_ind))), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxmed, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxmed, V_ind)),[],2)/sqrt(length(intersect(FIxmed, V_ind))), 'r');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Med visual: ' num2str(length(intersect(SIxmed, V_ind))) ' Hit; ' num2str(length(intersect(FIxmed, V_ind))) ' FA'])

subplot(2,3,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxlong, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxlong, V_ind)),[],2)/sqrt(length(intersect(SIxlong, V_ind))), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxlong, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxlong, V_ind)),[],2)/sqrt(length(intersect(FIxlong, V_ind))), 'r');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Long visual: ' num2str(length(intersect(SIxlong, V_ind))) ' Hit; ' num2str(length(intersect(FIxlong, V_ind))) ' FA'])

subplot(2,3,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxshort, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxshort, AV_ind)),[],2)/sqrt(length(intersect(SIxshort, AV_ind))), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxshort, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxshort, AV_ind)),[],2)/sqrt(length(intersect(FIxshort, AV_ind))), 'r');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Short auditory: ' num2str(length(intersect(SIxshort, AV_ind))) ' Hit; ' num2str(length(intersect(FIxshort, AV_ind))) ' FA'])

subplot(2,3,5)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxmed, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxmed, AV_ind)),[],2)/sqrt(length(intersect(SIxmed, AV_ind))), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxmed, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxmed, AV_ind)),[],2)/sqrt(length(intersect(FIxmed, AV_ind))), 'r');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Med auditory: ' num2str(length(intersect(SIxmed, AV_ind))) ' Hit; ' num2str(length(intersect(FIxmed, AV_ind))) ' FA'])

subplot(2,3,6)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxlong, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxlong, AV_ind)),[],2)/sqrt(length(intersect(SIxlong, AV_ind))), 'k');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxlong, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxlong, AV_ind)),[],2)/sqrt(length(intersect(FIxlong, AV_ind))), 'r');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Long auditory: ' num2str(length(intersect(SIxlong, AV_ind))) ' Hit; ' num2str(length(intersect(FIxlong, AV_ind))) ' FA'])
alignYaxes
suptitle([date ' ' mouse ' ' runstr '- align to previous stim on'])
print([fnout 'prevStimOn_align_SE_AVsep_SML.pdf'], '-dpdf')

figure;
tt = [-30:74].*(1000/30);
alignYaxes
subplot(2,3,1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxshort, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxshort, V_ind)),[],2)/sqrt(length(intersect(SIxshort, V_ind))), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxshort, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxshort, AV_ind)),[],2)/sqrt(length(intersect(SIxshort, AV_ind))), 'k');
xlim([-500 1500])
hold on
vline(0,'k')
title(['Short hit: ' num2str(length(intersect(SIxshort, V_ind))) ' Visual; ' num2str(length(intersect(FIxshort, AV_ind))) ' Auditory'])

subplot(2,3,2)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxmed, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxmed, V_ind)),[],2)/sqrt(length(intersect(SIxmed, V_ind))), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxmed, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxmed, AV_ind)),[],2)/sqrt(length(intersect(SIxmed, AV_ind))), 'k');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Med hit: ' num2str(length(intersect(SIxmed, V_ind))) ' Visual; ' num2str(length(intersect(FIxmed, AV_ind))) ' Auditory'])

subplot(2,3,3)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxlong, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxlong, V_ind)),[],2)/sqrt(length(intersect(SIxlong, V_ind))), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(SIxlong, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(SIxlong, AV_ind)),[],2)/sqrt(length(intersect(SIxlong, AV_ind))), 'k');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Long hit: ' num2str(length(intersect(SIxlong, V_ind))) ' Visual; ' num2str(length(intersect(SIxlong, AV_ind))) ' Auditory'])

subplot(2,3,4)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxshort, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxshort, V_ind)),[],2)/sqrt(length(intersect(FIxshort, V_ind))), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxshort, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxshort, AV_ind)),[],2)/sqrt(length(intersect(FIxshort, AV_ind))), 'k');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Short FA: ' num2str(length(intersect(FIxshort, V_ind))) ' Visual; ' num2str(length(intersect(FIxshort, AV_ind))) ' Auditory'])

subplot(2,3,5)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxmed, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxmed, V_ind)),[],2)/sqrt(length(intersect(FIxmed, V_ind))), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxmed, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxmed, AV_ind)),[],2)/sqrt(length(intersect(FIxmed, AV_ind))), 'k');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Med FA: ' num2str(length(intersect(FIxmed, V_ind))) ' Visual; ' num2str(length(intersect(FIxmed, AV_ind))) ' Auditory'])

subplot(2,3,6)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxlong, V_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxlong, V_ind)),[],2)/sqrt(length(intersect(FIxlong, V_ind))), 'g');
hold on
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,intersect(FIxlong, AV_ind)),2), nanstd(DataDFoverFavg(:,intersect(FIxlong, AV_ind)),[],2)/sqrt(length(intersect(FIxlong, AV_ind))), 'k');
hold on
xlim([-500 1500])
hold on
vline(0,'k')
title(['Long FA: ' num2str(length(intersect(FIxlong, V_ind))) ' Visual; ' num2str(length(intersect(FIxlong, AV_ind))) ' Auditory'])
suptitle([date ' ' mouse ' ' runstr '- align to previous stim on'])
print([fnout 'prevStimOn_align_SEsep_AV_SML.pdf'], '-dpdf')

%% Align data to previous stim on
Data = zeros(105,size(dataTC,2),ntrials);
DataDF = zeros(105,size(dataTC,2),ntrials);
DataDFoverF = zeros(105,size(dataTC,2),ntrials);
if cTargetOn(end,1)+30 > size(dataTC,1)
    for itrial = 1:ntrials-1
        if strcmp(input.trialOutcomeCell(itrial), 'failure')
            Data(:,:,itrial) = dataTC(cStimOn(itrial)-30:cStimOn(itrial)+74,:);
        else
            Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+74,:);
        end
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
else
    for itrial = 1:ntrials
        if strcmp(input.trialOutcomeCell(itrial), 'failure')
            Data(:,:,itrial) = dataTC(cStimOn(itrial)-30:cStimOn(itrial)+74,:);
        else
            Data(:,:,itrial) = dataTC(cTargetOn(itrial)-30:cTargetOn(itrial)+74,:);
        end
        DataDF(:,:,itrial) = bsxfun(@minus, Data(:,:,itrial), mean(Data(1:30,:,itrial),1));
        DataDFoverF(:,:,itrial) = bsxfun(@rdivide, DataDF(:,:,itrial), mean(Data(1:30,:,itrial),1));
    end
end

DataDFoverFavg = squeeze(mean(DataDFoverF,2));
tt = [-30:74].*(1000/30);
figure;
FIx = find(strcmp(input.trialOutcomeCell, 'failure'));
SIx = find(strcmp(input.trialOutcomeCell, 'success'));
MIx = find(strcmp(input.trialOutcomeCell, 'ignore'));
n = ceil(sqrt(length(cycles)));
if rem(n.^2, length(cycles))>n
    n2 = n-1;
else
    n2 = n;
end
for icyc = 1:length(cycles)
    FIxcyc = intersect(find(tCyclesOn==icyc), FIx);
    MIxcyc = intersect(find(tCyclesOn==icyc), MIx);
    RIxcyc = setdiff(find(tCyclesOn>icyc),MIxcyc);
    SIxcyc = intersect(find(tCyclesOn==icyc),SIx);
    Fb1Ix = intersect(V_ind, FIxcyc);
    Fb2Ix = intersect(AV_ind, FIxcyc);
    Rb1Ix = intersect(V_ind, RIxcyc);
    Rb2Ix = intersect(AV_ind, RIxcyc);
    Mb1Ix = intersect(V_ind, MIxcyc);
    Mb2Ix = intersect(AV_ind, MIxcyc);
    Sb1Ix = intersect(V_ind, SIxcyc);
    Sb2Ix = intersect(AV_ind, SIxcyc);
    Datacyc = zeros(105,size(dataTC,2),ntrials);
    DataDFcyc = zeros(105,size(dataTC,2),ntrials);
    DataDFoverFcyc = zeros(105,size(dataTC,2),ntrials);
    for itrial = 1:ntrials
        Datacyc(:,:,itrial) = dataTC(cLeverDown(itrial)+((icyc-1)*11)-30:cLeverDown(itrial)+((icyc-1)*11)+74,:);
        DataDFcyc(:,:,itrial) = bsxfun(@minus, Datacyc(:,:,itrial), mean(Datacyc(1:30,:,itrial),1));
        DataDFoverFcyc(:,:,itrial) = bsxfun(@rdivide, DataDFcyc(:,:,itrial), mean(Datacyc(1:30,:,itrial),1));
    end
    DataDFoverFcycavg = squeeze(mean(DataDFoverFcyc,2));
    if icyc == 1
        Hvc = figure;
        suptitle('Visual trials: False Alarm vs Correct reject')
        Hac = figure;
        suptitle('Auditory trials: False Alarm vs Correct reject')
        Hvs = figure;
        suptitle('Visual trials: False Alarm vs Success')
        Has = figure;
        suptitle('Auditory trials: False Alarm vs Success')
    end
    figure(Hvc)    
    subplot(n,n2,icyc)
    if length(Fb1Ix)>2
        shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
    else
        plot(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), 'r');
    end
    hold on
    if length(Rb1Ix)>2
        shadedErrorBar(tt, nanmean(DataDFoverFcycavg(:,Rb1Ix),2), nanstd(DataDFoverFcycavg(:,Rb1Ix),[],2)/sqrt(length(Rb1Ix)), 'k');
    else
        plot(tt, nanmean(DataDFoverFcycavg(:,Rb1Ix),2),'k');
    end
    xlim([-1000 1500])
    ylim([-.02 .05])
    hold on
    vline(0,'k')
    title([num2str(icyc) ' cyc:'  num2str(length(Rb1Ix)) ' CR; ' num2str(length(Fb1Ix)) ' FA']) 
    
    figure(Hac)    
    subplot(n,n2,icyc)
    if length(Fb2Ix)>2
        shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
    else
        plot(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), 'r');
    end
    hold on
    if length(Rb2Ix)>2
        shadedErrorBar(tt, nanmean(DataDFoverFcycavg(:,Rb2Ix),2), nanstd(DataDFoverFcycavg(:,Rb2Ix),[],2)/sqrt(length(Rb2Ix)), 'k');
    else
        plot(tt, nanmean(DataDFoverFcycavg(:,Rb2Ix),2),'k');
    end
    xlim([-1000 1500])
    ylim([-.02 .05])
    hold on
    vline(0,'k')
    title([num2str(icyc) ' cyc: ' num2str(length(Rb2Ix)) ' CR; ' num2str(length(Fb2Ix)) ' FA']) 
    
    figure(Hvs)    
    subplot(n,n2,icyc)
    if length(Fb1Ix)>2
        shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb1Ix),2), nanstd(DataDFoverFavg(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
    else
        plot(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), 'r');
    end
    hold on
    if length(Sb1Ix)>2
        shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
    else
        plot(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2),'k');
    end
    xlim([-1000 1500])
    ylim([-.02 .05])
    hold on
    vline(0,'k')
    title([num2str(icyc) ' cyc: ' num2str(length(Sb1Ix)) ' Hit; ' num2str(length(Fb2Ix)) ' FA']) 
    
    figure(Has)    
    subplot(n,n2,icyc)
    if length(Fb2Ix)>2
        shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), nanstd(DataDFoverFavg(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
    else
        plot(tt, nanmean(DataDFoverFavg(:,Fb2Ix),2), 'r');
    end
    hold on
    if length(Sb2Ix)>2
        shadedErrorBar(tt, nanmean(DataDFoverFcycavg(:,Sb2Ix),2), nanstd(DataDFoverFcycavg(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k');
    else
        plot(tt, nanmean(DataDFoverFcycavg(:,Sb2Ix),2),'k');
    end
    xlim([-1000 1500])
    ylim([-.02 .05])
    hold on
    vline(0,'k')
    title([num2str(icyc) ' cyc: ' num2str(length(Sb2Ix)) ' Hit; ' num2str(length(Fb2Ix)) ' FA']) 
end
figure(Hvc)
print([fnout 'prevStimOn_align_CF_Visual_bycyc.pdf'], '-dpdf')
figure(Hac)
print([fnout 'prevStimOn_align_CF_Auditory_bycyc.pdf'], '-dpdf')
figure(Hvs)
print([fnout 'prevStimOn_align_SF_Visual_bycyc.pdf'], '-dpdf')
figure(Has)
print([fnout 'prevStimOn_align_SF_Auditory_bycyc.pdf'], '-dpdf')


%% Align data to target on
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
n = length(Dirs);
for idir = 2:n
    subplot(ceil(n/2),2,idir-1)
    ind1 = intersect(Sb1Ix, find(tGratingDirectionDeg==Dirs(idir)));
    ind2 = intersect(Mb1Ix, find(tGratingDirectionDeg==Dirs(idir)));
    if length(ind1)>2
        shadedErrorBar(tt, nanmean(DataDFoverFavg(:,ind1),2), nanstd(DataDFoverFavg(:,ind1),[],2)./length(ind1), '-k');
        hold on;
    else
        plot(tt, nanmean(DataDFoverFavg(:,ind1),2), '-k')
        hold on;
    end
    if length(ind2)>2
        shadedErrorBar(tt, nanmean(DataDFoverFavg(:,ind2),2), nanstd(DataDFoverFavg(:,ind2),[],2)./length(ind2), '-r');
        hold on;
    else
        plot(tt, nanmean(DataDFoverFavg(:,ind2),2), '-r')
        hold on;
    end
    ylim([-0.02 0.04])
    xlim([-500 1500])
    vline([0 550],'--k')
    hold on
    title([num2str(Dirs(idir)) 'deg- success: ' num2str(length(ind1)) '; miss: ' num2str(length(ind2))])
end
subplot(ceil(n/2),2,idir)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb1Ix),2), nanstd(DataDFoverFavg(:,Sb1Ix),[],2)./length(Sb1Ix), '-k');
hold on;
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Mb1Ix),2), nanstd(DataDFoverFavg(:,Mb1Ix),[],2)./length(Mb1Ix), '-r');
xlim([-500 1500])
ylim([-0.02 0.04])
vline([0 550],'--k')
hold on
title('All visual')
subplot(ceil(n/2),2,idir+1)
shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Sb2Ix),2), nanstd(DataDFoverFavg(:,Sb2Ix),[],2)./length(Sb2Ix), '-k');
hold on;
if length(Mb2Ix) > 2
    shadedErrorBar(tt, nanmean(DataDFoverFavg(:,Mb2Ix),2), nanstd(DataDFoverFavg(:,Mb2Ix),[],2)./length(Mb2Ix), '-r');
else
    plot(tt, nanmean(DataDFoverFavg(:,Mb2Ix),2), '-r')
end
xlim([-500 1500])
ylim([-0.02 0.04])
vline([0 550],'--k')
hold on
title('All auditory')
suptitle([date ' ' mouse ' ' runstr '- align to target'])
print([fnout 'target_align_SM_AV.pdf'], '-dpdf')

reactIx = {};
reactIx{1,1} = intersect(V_ind, intersect(find(reactTime>100), find(reactTime<250)));
reactIx{2,1} = intersect(V_ind, intersect(find(reactTime>300), find(reactTime<400)));
figure;
colmat = strvcat('m', 'b');
nbin = 2;
subplot(2,1,1)
y = -.01:.01:.03;
for ibin = 1:2
    shadedErrorBar(tt, nanmean(DataDFoverFavg(:,reactIx{ibin,1}),2), nanstd(DataDFoverFavg(:,reactIx{ibin,1}),[],2)./length(reactIx{ibin,1}), colmat(ibin,:));
    hold on
    vline(mean(reactTime(1,reactIx{ibin,1}),2),colmat(ibin,:));
    hold on
end
vline([-660; -330; 0; 330],'--k')
xlim([-1000 1500])
title(['Visual: ' num2str(length(reactIx{1,1})) ' ' num2str(length(reactIx{2,1})) ' trials'])

subplot(2,1,2)
reactIx{1,2} = intersect(AV_ind, intersect(find(reactTime>100), find(reactTime<200)));
reactIx{2,2} = intersect(AV_ind, intersect(find(reactTime>300), find(reactTime<400)));
for ibin = 1:2
    shadedErrorBar(tt, nanmean(DataDFoverFavg(:,reactIx{ibin,2}),2), nanstd(DataDFoverFavg(:,reactIx{ibin,2}),[],2)./length(reactIx{ibin,2}), colmat(ibin,:));
    hold on
    vline(mean(reactTime(1,reactIx{ibin,2}),2),colmat(ibin,:));
    hold on
end
vline([-660; -330; 0; 330],'--k')
xlim([-1000 1500])
title(['Auditory: ' num2str(length(reactIx{1,2})) ' ' num2str(length(reactIx{2,2})) ' trials']) 
suptitle([date ' ' mouse ' ' runstr '- fast (magenta) and slow(blue) react times- all success'])
print([fnout 'target_align_byreact.pdf'], '-dpdf')
