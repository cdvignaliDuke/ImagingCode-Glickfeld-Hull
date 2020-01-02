%% expt info
fnout = 'Z:\home\ashley\Analysis\two-photon tinkering\scanboxAlign2Mworks';
rc = behavConstsAV;
subnum = '999';
mouse = 'Other\2P testing';
expDate = '190818';
runs = {'002','003'};
expTimes = {'1506','1510'};
nexp = length(runs);
exptType = {'flashingstim2Pframes','visstimret'};
%% get data
mw = cell(1,nexp);
data = cell(1,nexp);
for irun = 1:nexp
    runFolder = runs{irun};
    runTime = expTimes{irun};
    fname = [runFolder '_000_000'];
    fn = fullfile(rc.ashleyData,mouse, expDate,runFolder);
    cd(fn)
    mw{irun} = loadMworksFile(subnum,expDate,runTime);
    load([fname '.mat'])
    data{irun} = squeeze(sbxread(fname,0,info.config.frames));
end
rows = cellfun(@(x) x(101:end,:,:), data, 'unif',0);
tc = cellfun(@(x) squeeze(mean(mean(x,1),2)), rows, 'unif',0);

%% vis stim ret align to stim on
visStimRetInd = 2;

ip = mw{visStimRetInd};
d = tc{visStimRetInd};

on = ip.nScansOn;
off = ip.nScansOff;

nfr = length(d);
ntrials = floor(nfr./(on+off));
d = d(1:ntrials*(on+off));

d_tr = reshape(d,on+off,ntrials);

tFr = -off+1:on;

setFigParams4Print('portrait');figure;
suptitle('visStimRet')
subplot 211
plot(tFr,d_tr)
hold on
plot(tFr, mean(d_tr,2),'k-','LineWidth',2)
figXAxis([],'frames from onset',[])
figYAxis([],'F',[])
subplot 212
plot(tFr,d_tr)
hold on
plot(tFr, mean(d_tr,2),'k-','LineWidth',2)
figXAxis([],'frames from onset',[-1 5])
figYAxis([],'F',[330 340])
title('zoomed in')

print(fullfile(fnout,[expDate '_visStimRet']),'-dpdf','-fillpage')

%% stim test align to stim on
ip = mw{5};
d = tc{5};

on = ip.nScansOn;
off = ip.nScansOff;

nfr = length(d);
ntrials = floor(nfr./(on+off));
d = d(1:ntrials*(on+off));

d_tr = reshape(d,on+off,ntrials);

tFr = -off+1:on;

setFigParams4Print('portrait');figure;
suptitle('stimTest')
subplot 211
plot(tFr,d_tr)
figXAxis([],'frames from onset',[])
figYAxis([],'F',[])
subplot 212
plot(tFr,d_tr)
figXAxis([],'frames from onset',[-1 5])
figYAxis([],'F',[975 982])
title('zoomed in')

print(fullfile(fnout,'stimTest'),'-dpdf','-fillpage')

%% flashingStim_2P_Frames align to first stim and catch stim
fsInd = 1;

preframes = 30;
postframes = 60;
tt = (1:(preframes+postframes))-preframes;

ip = mw{fsInd};
d = tc{fsInd};

trStart = cell2mat(ip.cFirstStim);
ntrials = length(trStart)-1;
trStart = trStart(1:ntrials);
d_tr = nan(preframes+postframes,ntrials);
for itrial = 1:ntrials
   ind = trStart(itrial);
   indfr = ind-preframes+1:ind+postframes;
   d_tr(:,itrial) = d(indfr);    
end

targetStart = cell2mat(ip.cTargetOn);
ntrials = length(targetStart)-1;
targetStart = targetStart(1:ntrials);
d_target = nan(preframes+postframes,ntrials);
for itrial = 1:ntrials
   ind = targetStart(itrial);
   indfr = ind-preframes+1:ind+postframes;
   d_target(:,itrial) = d(indfr);    
end

catchStart = celleqel2mat_padded(ip.cCatchOn);
catchStart = catchStart(~isnan(catchStart));
ntrials = length(catchStart);
d_catch = nan(preframes+postframes,ntrials);
for itrial = 1:ntrials
   ind = catchStart(itrial);
   indfr = ind-preframes+1:ind+postframes;
   d_catch(:,itrial) = d(indfr);    
end

setFigParams4Print('landscape');figure;
suptitle('FS 2P Frames')
subplot 221
plot(tt,d_target)
hold on
h = plot(tt,mean(d_target,2),'k');
h.LineWidth = 1;
figXAxis([],'frames from cTargetStim',[])
figYAxis([],'F',[])
subplot 222
plot(tt,d_catch)
hold on
h = plot(tt,mean(d_catch,2),'k');
h.LineWidth = 1;
figXAxis([],'frames from cCatchOn',[])
figYAxis([],'F',[])
subplot 223
plot(tt,d_tr)
hold on
h = plot(tt,mean(d_tr,2),'k');
h.LineWidth = 1;
figXAxis([],'frames from cFirstStim',[])
figYAxis([],'F',[])
subplot 224
h = plot(tt,mean(d_target,2),'k');
hold on
h = plot(tt,mean(d_catch,2),'b');
h = plot(tt,mean(d_tr,2),'r');
h.LineWidth = 1;
figXAxis([],'frames from stim on',[-10 10])
figYAxis([],'F',[])
legend({'Target Stim';'Catch Stim';'First Stim'},'location','northwest')

print(fullfile(fnout,[expDate '_FS2PFrames_withCatch']),'-dpdf','-fillpage')


%% holdanddetect_2P_Frames align to target
preframes = 50;
postframes = 50;
ip = mw{3};
d = tc{3};

trStart = cell2mat(ip.cTargetOn);
ntrials = length(trStart)-1;
trStart = trStart(1:ntrials);

d_tr = nan(preframes+postframes,ntrials);
for itrial = 1:ntrials
   ind = trStart(itrial);
   indfr = ind-preframes+1:ind+postframes;
   d_tr(:,itrial) = d(indfr);    
end

setFigParams4Print('portrait');figure;
suptitle('HD 2P Frames')
subplot 211
plot(tFr,d_tr)
figXAxis([],'frames from onset',[])
figYAxis([],'F',[])
subplot 212
plot(tFr,d_tr)
figXAxis([],'frames from onset',[-1 5])
figYAxis([],'F',[970 982])
title('zoomed in')

print(fullfile(fnout,'HD2PFrames_fast'),'-dpdf','-fillpage')

%% holdanddetect_2P_Frames align to target
preframes = 50;
postframes = 50;
ip = mw{4};
d = tc{4};

trStart = cell2mat(ip.cTargetOn);
ntrials = length(trStart)-1;
trStart = trStart(1:ntrials);

d_tr = nan(preframes+postframes,ntrials);
for itrial = 1:ntrials
   ind = trStart(itrial);
   indfr = ind-preframes+1:ind+postframes;
   d_tr(:,itrial) = d(indfr);    
end

setFigParams4Print('portrait');figure;
suptitle('HD 2P Frames')
subplot 211
plot(tFr,d_tr)
figXAxis([],'frames from onset',[])
figYAxis([],'F',[])
subplot 212
plot(tFr,d_tr)
figXAxis([],'frames from onset',[-1 5])
figYAxis([],'F',[970 982])
title('zoomed in')

print(fullfile(fnout,'HD2PFrames_long'),'-dpdf','-fillpage')
