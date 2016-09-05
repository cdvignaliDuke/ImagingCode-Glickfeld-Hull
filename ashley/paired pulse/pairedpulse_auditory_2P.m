%experiment info
pairedPulseAud_datasets
iexp = 2;
mouse = expt(iexp).mouse;
SubNum = expt(iexp).SubNum;
expdate = expt(iexp).date;
mw_time = expt(iexp).time_mat;
folder = expt(iexp).folder;
exprun = expt(iexp).runs;

frameRate = expt(iexp).frame_rate;

fnin = fullfile('Z:\Analysis',mouse,folder,expdate,exprun);
fnout = fullfile('Z:\Analysis',mouse,folder,expdate);
mkdir(fnout);

%load data
load(fullfile(fnin,'Timecourses.mat'));
load(['Y:\home\andrew\Behavior\Data\data-i' SubNum '-' expdate '-' mw_time])

if iexp == 2
    data = data_TC;
else
    data = dataTimecourse.dataTCsub;
end

%% set params for figures
set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);


%important variables
AVtrials = find(cell2mat(input.tBlock2TrialNumber) == 0);
Atrials = find(cell2mat(input.tBlock2TrialNumber) == 1);
ntrials = input.trialSinceReset;
ncells = size(data,2);
preTrialFrames = 30;
stimON = input.nFramesOn;
stimOFF = input.nFramesOff;
trialFrames = (stimON+stimOFF)*2;
pulseFrames = stimON+stimOFF;
trialStart_ind = cell2mat(input.cLeverDown);
pulseStart_ind = zeros(length(trialStart_ind)*2,1);
pulseStart_ind(1:2:length(trialStart_ind)*2) = trialStart_ind;
pulseStart_ind(2:2:length(trialStart_ind)*2) = trialStart_ind+(pulseFrames);
pulse1_win = preTrialFrames+1:preTrialFrames+stimON;
pulse2_win = preTrialFrames+pulseFrames:preTrialFrames+(2*stimON)+stimOFF;

isblk2 = cell2mat(input.tBlock2TrialNumber);
pulseType = zeros(ntrials*2,1);
pulseType(1:2:ntrials*2) = isblk2;
pulseType(2:2:ntrials*2) = isblk2+2;
AV_ind = find(pulseType == 0);
Aonly_ind = find(pulseType == 1);
Vonly_ind = find(pulseType == 2);
blank_ind = find(pulseType == 3);

%find trial-by-trial dF/F
dFoverF = zeros(preTrialFrames+trialFrames,ncells,ntrials);
for i = 1:ntrials
    F = nanmean(data(trialStart_ind(i)-preTrialFrames:trialStart_ind(i),:),1);
    dFoverF(:,:,i) = bsxfun(@rdivide, bsxfun(@minus, data(trialStart_ind(i)-preTrialFrames:trialStart_ind(i)+trialFrames-1,:),F),F);
end

%find pulse-by-pulse dF/F
dFoverF_pulse = zeros(preTrialFrames+pulseFrames,ncells,ntrials*2);
for i = 1:ntrials*2
    F = mean(data(pulseStart_ind(i)-preTrialFrames:pulseStart_ind(i),:),1);
    dFoverF_pulse(:,:,i) = bsxfun(@rdivide, bsxfun(@minus, data(pulseStart_ind(i)-preTrialFrames:pulseStart_ind(i)+pulseFrames-1,:),F),F);
end


%timecourses
ttFrames = -preTrialFrames:trialFrames-1;
ttMs = (ttFrames/frameRate)*1000;
AV_tc = mean(mean(dFoverF(:,:,AVtrials),3),2);
A_tc = mean(mean(dFoverF(:,:,Atrials),3),2);

AV_err = std(squeeze(mean(dFoverF(:,:,AVtrials),2)),[],2)./sqrt(length(AVtrials));
A_err = std(squeeze(mean(dFoverF(:,:,Atrials),2)),[],2)./sqrt(length(Atrials));

figure;
AV = shadedErrorBar(ttMs,AV_tc,AV_err,'c');
hold on
A = shadedErrorBar(ttMs,A_tc,A_err,'b');
hold on
vline(0,'k-')
hold on
vline(((pulseFrames)/frameRate)*1000,'k-')



%timecourses of indiv. pulses
ttMs_pulse = ((-preTrialFrames:pulseFrames-1)/frameRate)*1000;
AV_pulse_tc = mean(mean(dFoverF_pulse(:,:,AV_ind),3),2);
Aonly_pulse_tc = mean(mean(dFoverF_pulse(:,:,Aonly_ind),3),2);
Vonly_pulse_tc = mean(mean(dFoverF_pulse(:,:,Vonly_ind),3),2);
blank_pulse_tc = mean(mean(dFoverF_pulse(:,:,blank_ind),3),2);

AV_pulse_err = std(squeeze(mean(dFoverF_pulse(:,:,AV_ind),2)),[],2)./sqrt(length(AV_ind));
Aonly_pulse_err = std(squeeze(mean(dFoverF_pulse(:,:,Aonly_ind),2)),[],2)./sqrt(length(Aonly_ind));
Vonly_pulse_err = std(squeeze(mean(dFoverF_pulse(:,:,Vonly_ind),2)),[],2)./sqrt(length(Vonly_ind));
blank_pulse_err = std(squeeze(mean(dFoverF_pulse(:,:,blank_ind),2)),[],2)./sqrt(length(blank_ind));

figure;
subplot(1,2,1)
AV = shadedErrorBar(ttMs_pulse,AV_pulse_tc,AV_pulse_err,'c');
hold on
A = shadedErrorBar(ttMs_pulse,Aonly_pulse_tc,Aonly_pulse_err,'r');
hold on
V = shadedErrorBar(ttMs_pulse,Vonly_pulse_tc,Vonly_pulse_err,'b');
hold on
blank = shadedErrorBar(ttMs_pulse,blank_pulse_tc,blank_pulse_err,'m');
hold on
vline(0,'k--')
legend([AV.mainLine,A.mainLine,V.mainLine,blank.mainLine],{'AV';'A only';'V only';'blank'})
xlabel('time (ms)');
ylabel('dF/F');

subplot(1,2,2)
AV = plot(ttMs_pulse,AV_pulse_tc,'c');
hold on
A = plot(ttMs_pulse,Aonly_pulse_tc,'r');
hold on
V = plot(ttMs_pulse,Vonly_pulse_tc,'b');
hold on
blank = plot(ttMs_pulse,blank_pulse_tc,'m');
hold on
vline(0,'k--')
legend([AV,A,V,blank],{'AV';'A only';'V only';'blank'})
xlabel('time (ms)');
ylabel('dF/F');

%% plot individual cells' traces for each condition

figure
subplot(2,2,1)
plot(ttMs_pulse,squeeze(mean(dFoverF_pulse(:,:,AV_ind),3)))
hold on
title('visual + auditory')
subplot(2,2,2)
plot(ttMs_pulse,squeeze(mean(dFoverF_pulse(:,:,Vonly_ind),3)))
title('visual only')
hold on
subplot(2,2,3)
plot(ttMs_pulse,squeeze(mean(dFoverF_pulse(:,:,Aonly_ind),3)))
hold on
title('auditory only')
subplot(2,2,4)
plot(ttMs_pulse,squeeze(mean(dFoverF_pulse(:,:,blank_ind),3)))
hold on
title('blank')
for i = 1:4
    subplot(2,2,i)
    xlabel('time (ms)')
    ylabel('dF/F')
    ylim([-0.1 1])
    vline([0 500],'k--')
end
print([fnout '\' expt(iexp).img_loc{1} 'cellresp_bytrialtype'],'-dpdf');

%% find responses to trial start for each 
resp_win = 3; % number of frames before and after peak
prestim_win = preTrialFrames+(0.1*frameRate);
stim_win = prestim_win:prestim_win+(0.5*frameRate);
cells.name = {'AV';'V only';'A only';'blank'};
cells.tc = [{mean(dFoverF_pulse(:,:,AV_ind),3)},{mean(dFoverF_pulse(:,:,Vonly_ind),3)},{mean(dFoverF_pulse(:,:,Aonly_ind),3)},{mean(dFoverF_pulse(:,:,blank_ind),3)}];
[cells.peakval cells.peakind] = max(mean(dFoverF_pulse(stim_win,:,:),3),[],1);
cells.peakind = cells.peakind+prestim_win;

cells.respwin_ind = linspaceNDim(cells.peakind'-(resp_win-1),cells.peakind'+(resp_win-1),(resp_win*2)-1);

cells.resp = cellfun(@(x) mean(x(cells.respwin_ind,:),1),cells.tc,'unif',false);

cells.mi = bsxfun(@rdivide,bsxfun(@plus,cells.resp{2},cells.resp{3}),cells.resp{1});

bin_width = 0.25;
bin_edges = [min(cells.mi)-1 (-2:bin_width:2) max(cells.mi)+1];

[cells.bin_n cells.bin_ind] = histc(cells.mi,bin_edges);
cells.bin_avg = zeros(1,max(cells.bin_ind)+1);
for i = 1:max(cells.bin_ind)
    cells.bin_avg(i) = mean(cells.mi(find(cells.bin_ind == i)));
end

figure;
hstgrm = bar(cells.bin_avg(~isnan(cells.bin_avg)),cells.bin_n(~isnan(cells.bin_avg)));
hstgrm.Parent.XTick = -2:2;%chop(cells.bin_avg(~isnan(cells.bin_avg)),2);
text(-4,10,{['mean = ' num2str(mean(cells.mi))]; ['std = ' num2str(std(cells.mi))]});
hstgrm.Parent.YLim = [0 max(cells.bin_n)+1];
hstgrm.Parent.XLim = [-5 5];
hstgrm.Parent.XLabel.String = 'MI';
hstgrm.Parent.YLabel.String = 'n cells';
title({expt(iexp).img_loc{1}; ' - MI = (A+V)/AV'})
print([fnout '\' expt(iexp).img_loc{1} 'MI_AplusVoverAV'],'-dpdf');


%% find responsive cells
V = dFoverF_pulse(:,:,Vonly_ind);
A = dFoverF_pulse(:,:,Aonly_ind);

[h_V p_V] = ttest(mean(V(1:preTrialFrames,:,:),1),mean(V(stim_win,:,:),1),'dim',3,'tail','left');
[h_A p_A] = ttest(mean(A(1:preTrialFrames,:,:),1),mean(A(stim_win,:,:),1),'dim',3,'tail','left');

Vis_ind = find(h_V);
Aud_ind = find(h_A);
Both_ind = intersect(Vis_ind,Aud_ind);


figure;
subplot(2,2,1)
plot(ttMs_pulse,squeeze(mean(dFoverF_pulse(:,Vis_ind,Vonly_ind),3)))
title('visual only')
hold on
subplot(2,2,2)
plot(ttMs_pulse,squeeze(mean(dFoverF_pulse(:,Aud_ind,Aonly_ind),3)))
title('auditory only')
hold on
subplot(2,2,3)
plot(ttMs_pulse,squeeze(mean(dFoverF_pulse(:,setdiff(1:60,Vis_ind),Vonly_ind),3)))
title('visual only - non-resp')
hold on

%% scatter V vs AV
figure;
h = scatter(cells.resp{1}, cells.resp{2},'k')
hold on
plot([-1:.5:1],[-1:.5:1],'k--')
xlim([-.1 .2]);
ylim([-.1 .2]);
axis square
hold on
xlabel('AV')
ylabel('V')
title(expt(iexp).img_loc{1})
print([fnout '\' expt(iexp).img_loc{1} 'VvsAVscat'],'-dpdf');





