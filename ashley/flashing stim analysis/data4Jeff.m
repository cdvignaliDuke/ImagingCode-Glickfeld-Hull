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
iexp = 1;

%%
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
frameRate = 30;
stimOnS = 0.1;
stimOnFr = stimOnS*frameRate;

visStimDelayFr = 7; %trial start is cLeverDown, not cFirstStim
respwin = (nBaselineFr+visStimDelayFr):...
    (nBaselineFr+visStimDelayFr+stimOnFr-1);
basewin = (nBaselineFr-1):nBaselineFr+1;
respwin_tar = (nBaselineFr+visStimDelayFr-2):...
    (nBaselineFr+visStimDelayFr+stimOnFr-3);
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

firstStimTC = firstStimResp;
firstStimResp = nan(nBaselineFr*2,nCells,nTrials);
nAlignedResp = nan(nBaselineFr*2,nCells,nTrials);
for i = 1:nTrials
    ind = (((n(i)-1)*cycFr)-nBaselineFr):(((n(i)-1)*cycFr)+nBaselineFr-1);
    nAlignedResp(:,:,i) = r(ind,:,i);
    ind = 1:(2*nBaselineFr);
    firstStimResp(:,:,i) = r(ind,:,i);
end

firstStimTC = cat(3,firstStimTC,firstStimResp);
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
targetStimRespBL = targetStimResp - mean(targetStimResp(basewin_tar,:,:),1);
targetResp = squeeze(mean(targetStimRespBL(respwin_tar,:,:),1));

targetTC = targetStimRespBL;

r = dMiss.resp;
[~,nCells,nTrials] = size(r);

targetStimResp = nan(nBaselineFr*2,nCells,nTrials);
for i = 1:nTrials
    ind = 1:(2*nBaselineFr);
    targetStimResp(:,:,i) = r(ind,:,i);
end
targetStimRespBL = targetStimResp - mean(targetStimResp(basewin_tar,:,:),1);
targetResp = cat(2,targetResp,squeeze(mean(targetStimRespBL(respwin_tar,:,:),1)));
targetTC = cat(3,targetTC,targetStimRespBL);

%%
oris = unique(ori);
nori = length(oris);
df = (nori+1)-1;
respCellAlpha = 0.05./df;

firstResponsiveCells = logical(ttest(squeeze(mean(firstStimTC(respwin,:,:),1)),...
    squeeze(mean(firstStimTC(basewin,:,:),1)),...
    'dim',2,'tail','right','alpha',respCellAlpha));

targetResponsiveCells = nan(nCells,nori);
for i = 1:nori
    ind = ori == oris(i);
    r = squeeze(mean(targetTC(respwin_tar,:,ind),1));
    b = squeeze(mean(targetTC(basewin_tar,:,ind),1)); 
    targetResponsiveCells(:,i) = logical(ttest(r,b,...
    'dim',2,'tail','right','alpha',respCellAlpha));
end

allTargetResponsiveCells = logical(ttest(squeeze(mean(targetTC(respwin_tar,:,:),1)),...
    squeeze(mean(targetTC(basewin_tar,:,:),1)),...
    'dim',2,'tail','right','alpha',0.05));

isResponsiveCell = firstResponsiveCells | sum(targetResponsiveCells,2) > 0 | allTargetResponsiveCells;
%%
expt.firstBaseResp = firstResp;
expt.lastBaseResp = nAlignedRespBaselined;
expt.trialOutcome = trialOut;
expt.targetOrientation = ori;
expt.targetResp = targetResp;
expt.signifResponsiveCells = isResponsiveCell;
%%
fn = fullfile(rc.ashleyAnalysis,expt.mouse,'two-photon imaging',...
    expt.date);
if ~exist(fn,'dir')
    mkdir(fn)
end
save(fullfile(fn,'bxImgSum'),'expt')

leg = [];
figure
ind = firstResponsiveCells;
subplot 121
plot(mean(mean(firstStimTC(:,ind,:),3),2))
hold on
vline((0:cycFr:22)+nBaselineFr,'k:');
h = vline(nBaselineFr,'k:');
leg(1) = h;
vline(respwin(1),'r')
h = vline(respwin(end),'r');
leg(2) = h;
legend(leg,{'base stim on';'resp win'},'location','northwest')
figXAxis([],'Frame Number',[])
figYAxis([],'dF/F',[])
figAxForm
title(sprintf('First Stim Responsive Cells (%s)',num2str(sum(ind))))

ind = sum(targetResponsiveCells,2) > 0 | allTargetResponsiveCells;
subplot 122
plot(mean(mean(targetTC(:,ind,:),3),2))
hold on
vline((0:cycFr:22)+nBaselineFr,'k:');
h = vline(nBaselineFr,'k:');
vline(respwin_tar(1),'r')
h = vline(respwin_tar(end),'r');
figXAxis([],'Frame Number',[])
figYAxis([],'dF/F',[])
figAxForm
title(sprintf('Target Responsive Cells (%s)',num2str(sum(ind))))

print(fullfile(fn,'resp2stim'),'-dpdf','-fillpage')