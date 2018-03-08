clear all
close all
ds = 'FSAV_V1_100ms';
cellsOrDendrites = 1;
doCellLabel = 0;
%%
rc = behavConstsAV;
if strcmp(rc.name,'ashle') & strcmp(ds(1:3),'FSA')
    dataGroup = ds;
else
    dataGroup = [];
end
eval(dataGroup)

titleStr = ds;
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];

if cellsOrDendrites == 1
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary_cells' ds(5:end) '.mat']));    
elseif cellsOrDendrites == 2
    load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary_dendrites' ds(5:end) '.mat']));    
end        
if cellsOrDendrites == 1
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_startAlign_']); 
elseif cellsOrDendrites == 2
    fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_den_startAlign_']); 
end
   
%%
imouse = 2;
iexp = 2;
pressAlign = 1;
targetAlign = 2;
visualTrials = 1;
auditoryTrials = 2;
hitTrials = 1;
missTrials = 2;

%%
d = mouse(imouse).expt(iexp);

nBaselineFr = d.pre_event_frames;
cycFr = d.info.cyc_time;
nStimFr = nBaselineFr;
frameRate = expt(1).frame_rate;
stimOnS = 0.1;
stimOnFr = stimOnS*frameRate;

visStimDelayFr = 5; %trial start is cLeverDown, not cFirstStim
respwin = (nBaselineFr+visStimDelayFr):...
    (nBaselineFr+visStimDelayFr+stimOnFr-1);
basewin = (nBaselineFr-1):nBaselineFr+1;
respwin_tar = (nBaselineFr+visStimDelayFr-1):...
    (nBaselineFr+visStimDelayFr+stimOnFr-2);
basewin_tar = (nBaselineFr-1):nBaselineFr+1;
%%
expt = struct;
expt.mouse = d.mouse_name;
expt.date = d.date;

expt.info.respdata = 'cells x trials';
%%
dHit = d.align(pressAlign).av(visualTrials).outcome(hitTrials);
dMiss = d.align(pressAlign).av(visualTrials).outcome(missTrials);

r = dHit.resp;
n = dHit.tcyc;
nt = length(n);
[~,nCells,nTrials] = size(r);

firstStimResp = nan(nBaselineFr*2,nCells,nTrials);
nAlignedResp = nan(nBaselineFr*2,nCells,nTrials);
for i = 1:nTrials
    ind = (((n(i)-1)*cycFr)-nBaselineFr):(((n(i)-1)*cycFr)+nBaselineFr-1);
    nAlignedResp(:,:,i) = r(ind,:,i);
    ind = 1:(2*nBaselineFr);
    firstStimResp(:,:,i) = r(ind,:,i);
end

trialOut = repmat({'h'},1,nt);
firstResp = squeeze(mean(firstStimResp(respwin,:,:),1));
nAlignedRespBaselined = squeeze(mean(nAlignedResp(respwin,:,:),1)...
    - mean(nAlignedResp(basewin,:,:),1));

r = dMiss.resp;
n = dMiss.tcyc;
nt = length(n);
[~,nCells,nTrials] = size(r);

firstStimResp = nan(nBaselineFr*2,nCells,nTrials);
nAlignedResp = nan(nBaselineFr*2,nCells,nTrials);
for i = 1:nTrials
    ind = (((n(i)-1)*cycFr)-nBaselineFr):(((n(i)-1)*cycFr)+nBaselineFr-1);
    nAlignedResp(:,:,i) = r(ind,:,i);
    ind = 1:(2*nBaselineFr);
    firstStimResp(:,:,i) = r(ind,:,i);
end

ori = cat(2,dHit.targetOri,dMiss.targetOri);
trialOut = cat(2,trialOut,repmat({'m'},1,nt));
firstResp = cat(2,firstResp,squeeze(mean(firstStimResp(respwin,:,:),1)));
nAlignedRespBaselined = cat(2,nAlignedRespBaselined,...
    squeeze(mean(nAlignedResp(respwin,:,:),1)...
    - mean(nAlignedResp(basewin,:,:),1)));
%%
dHit = d.align(targetAlign).av(visualTrials).outcome(hitTrials);
dMiss = d.align(targetAlign).av(visualTrials).outcome(missTrials);

r = dHit.resp;
[~,nCells,nTrials] = size(r);

targetStimResp = nan(nBaselineFr*2,nCells,nTrials);
for i = 1:nTrials
    ind = 1:(2*nBaselineFr);
    targetStimResp(:,:,i) = r(ind,:,i);
end
targetResp = squeeze(mean(targetStimResp(respwin_tar,:,:),1));

r = dMiss.resp;
[~,nCells,nTrials] = size(r);

targetStimResp = nan(nBaselineFr*2,nCells,nTrials);
for i = 1:nTrials
    ind = 1:(2*nBaselineFr);
    targetStimResp(:,:,i) = r(ind,:,i);
end
targetResp = cat(2,targetResp,squeeze(mean(targetStimResp(respwin_tar,:,:),1)));

%%
expt.firstBaseResp = firstResp;
expt.lastBaseResp = nAlignedRespBaselined;
expt.trialOutcome = trialOut;
expt.targetOrientation = ori;
expt.targetResp = targetResp;
%%
fn = fullfile(rc.ashleyAnalysis,expt.mouse,'two-photon imaging',...
    expt.date);
mkdir(fn)
save(fullfile(fn,'bxImgSum'),'expt')
