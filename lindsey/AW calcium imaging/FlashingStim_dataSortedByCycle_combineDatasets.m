%combine two datasets to have 3 trial types - vis only, aud only, and
%vis+aud
runstr = runs(1,:);
if nrun>1
    for irun = 2:nrun
        runstr = [runstr '-' runs(irun,:)];
    end
end
fnout = fullfile('Z:\home\lindsey\Analysis\2P', mouse, date, [date '_' mouse '_' runstr '_']);
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
suptitle([date ' ' mouse ' ' runstr '- align to lever release'])
print([fnout 'release_align_SE_AV.pdf'], '-dpdf')

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
