%% expt info
fnout = 'Z:\Analysis\two-photon tinkering\scanboxAlign2Mworks';
rc = behavConstsAV;
subnum = '999';
mouse = 'Other\2P testing';
expDate = '170710';
runs = {'001';'002';'003';'004';'005'};
expTimes = {'1628';'1631';'1636';'1638';'1643'};
nexp = length(runs);
%% get data
mw = cell(1,nexp);
data = cell(1,nexp);
for irun = 1:nexp
    runFolder = runs{irun};
    runTime = expTimes{irun};
    fname = [runFolder '_000_000'];
    fn = fullfile(rc.ashleyAnalysis,mouse, expDate,runFolder);
    [mw{irun},data{irun}] = Load_SBXdataPlusMWorksData(subnum,expDate,runTime,mouse,runFolder,fname);
end
rows = cellfun(@(x) x(101:end,:,:), data, 'unif',0);
tc = cellfun(@(x) squeeze(mean(mean(x,1),2)), rows, 'unif',0);

%% vis stim ret align to stim on
ip = mw{1};
d = tc{1};

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
figXAxis([],'frames from onset',[])
figYAxis([],'F',[])
subplot 212
plot(tFr,d_tr)
figXAxis([],'frames from onset',[-1 5])
figYAxis([],'F',[970 975])
title('zoomed in')

print(fullfile(fnout,'visStimRet'),'-dpdf','-fillpage')


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

%% flashingStim_2P_Frames align to trial start stim on
preframes = 50;
postframes = 50;
ip = mw{2};
d = tc{2};

trStart = cell2mat(ip.cStimOn);
ntrials = length(trStart)-1;
trStart = trStart(1:ntrials);

d_tr = nan(preframes+postframes,ntrials);
for itrial = 1:ntrials
   ind = trStart(itrial);
   indfr = ind-preframes+1:ind+postframes;
   d_tr(:,itrial) = d(indfr);    
end

setFigParams4Print('portrait');figure;
suptitle('FS 2P Frames')
subplot 211
plot(tFr,d_tr)
figXAxis([],'frames from onset',[])
figYAxis([],'F',[])
subplot 212
plot(tFr,d_tr)
figXAxis([],'frames from onset',[-1 5])
figYAxis([],'F',[970 982])
title('zoomed in')

print(fullfile(fnout,'FS2PFrames'),'-dpdf','-fillpage')


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
