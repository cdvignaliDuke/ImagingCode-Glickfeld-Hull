ds = 'audDelay_V1_SOM';
eval(ds)
rc = behavConstsAV;

fn = fullfile(rc.ashleyAnalysis,'Expt Summaries',ds);
load(fullfile(fn,'dataSummary'))

doLoadPreviousTuning = false;

basewinS = 0.1;
respwinS = 0.1;

minVisOnlyMs = 8000;
nBoot = 10;

orisSmooth = 0:180;
tuningCenterFit = 90;
fitReliabilityCutoff = 22.5;

nBaselineFr = params.nBaselineMs.*params.frameRate./1000;
basewin = nBaselineFr-round(basewinS*params.frameRate):nBaselineFr;
respwin = (nBaselineFr+params.nFramesVisDelay):...
    (nBaselineFr+params.nFramesVisDelay+(respwinS.*params.frameRate));


im = 2;
iexp = 1;

tOris = ms(im).expt(iexp).tOrientation;
orientations = unique(tOris);
oriRespAlpha = 0.05./(length(orientations));
tuningCenterResp = find(orientations == tuningCenterFit);
nOri = length(orientations);
tDelay = ms(im).expt(iexp).tAudVisDelay;

visTC = ms(im).expt(iexp).tVisAlignResp;
audTC = ms(im).expt(iexp).tAudAlignResp;

visResp = squeeze(mean(visTC(respwin,:,:),1));
blResp = squeeze(mean(visTC(basewin,:,:),1));

[trialLength,nCells,~] = size(visTC);

visTuning = nan(nCells,nOri);
visTuningErr = nan(nCells,nOri);
nTrialsPerOri = zeros(1,nOri);
tVisRespOri = cell(1,nOri);
tVisBLOri = cell(1,nOri);
for i = 1:nOri
    ind = tOris == orientations(i) & tDelay >= 6000;
    visTuning(:,i) = mean(visResp(:,ind),2);
    visTuningErr(:,i) = ste(visResp(:,ind),2);
    nTrialsPerOri(i) = sum(ind);
    
    tVisRespOri{i} = visResp(:,ind);
    tVisBLOri{i} = blResp(:,ind);
end

oriAlpha = 0.05/nOri;
oriTTest = sum(cell2mat(cellfun(@(x,y) ...
    ttest(x,y,'dim',2,'tail','right','alpha',oriAlpha),...
    tVisRespOri,tVisBLOri,'unif',0)),2) > 0;
allTTest = ttest(cell2mat(tVisRespOri),cell2mat(tVisBLOri),...
    'dim',2,'tail','right');
isResponsive = oriTTest | allTTest;

if doLoadPreviousTuning
    load(fullfile(fn,'oriTuning.mat'))
else
    visTuning_resamp = nan(nCells,nOri,nBoot);
    for iboot = 1:nBoot
        for iori = 1:nOri
            ind = find(tOris == orientations(i) & tDelay >= 6000);   
            randInd = randsample(ind,length(ind),1);
            visTuning_resamp(:,iori,iboot) = mean(visResp(:,randInd),2);
        end
    end   
    [vonMisesFit_boot,~,fitReliability,R_square] = vonmisesReliableFit(...
        visTuning,visTuning_resamp,orientations,nBoot);
    vonMisesFit = squeeze(vonMisesFit_boot(:,1,:));
    R_square = R_square(1,:);
    save(fullfile(fn,'oriTuning'),'vonMisesFit','fitReliability','R_square')
end

%%
% isTuned = fitReliability > fitReliabilityCutoff;
isTuned = R_square > 0.9;

delays = unique(tDelay);
nDelay = length(delays);

visDelayTuning = cell(1,nDelay);
visDelayTuningErr = cell(1,nDelay);
visDelayTuningTTest = cell(1,nDelay);
visDelayRespTTest = cell(1,nDelay);
visTCTuning = cell(1,nDelay); 
visTCTuningErr = cell(1,nDelay);
nTrialsEaDelay = zeros(nOri,nDelay);
tVisTCTuning = cell(1,nDelay);
for i = 1:nDelay
    tuning = nan(nCells,nOri);
    tuningErr = nan(nCells,nOri);
    tc = nan(trialLength,nCells,nOri);
    tcErr = nan(trialLength,nCells,nOri);
    for iori = 1:nOri
        ind = tOris == orientations(iori) & tDelay == delays(i);
        tuning(:,iori) = mean(visResp(:,ind),2);
        tuningErr(:,iori) = ste(visResp(:,ind),2);
        tuning_ttest(:,iori) = ttest(...
            visResp(:,ind),blResp(:,ind),...
            'dim',2,'tail','right','alpha',oriRespAlpha);
        
        tc(:,:,iori) = mean(visTC(:,:,ind),3);
        tcErr(:,:,iori) = ste(visTC(:,:,ind),3);
        nTrialsEaDelay(iori,i) = sum(ind);
        tVisTCTuning{i}{iori} = visTC(:,:,ind);
    end    
    ind = tDelay == delays(i);
    alltuning_ttest = ttest(...
            visResp(:,ind),blResp(:,ind),...
            'dim',2,'tail','right','alpha',0.05);
    visDelayTuning{i} = tuning;
    visDelayTuningErr{i} = tuningErr;
    visDelayTuningTTest{i} = tuning_ttest;
    visDelayRespTTest{i} = alltuning_ttest;
    visTCTuning{i} = tc;
    visTCTuningErr{i} = tcErr;
end

vonMisesFits_delays = cell(1,nDelay);
for i = 1:nDelay
    tuning = visDelayTuning{i};
    y_fit = nan(length(orisSmooth),nCells);
    for iCell = 1:nCells
        [b, k, R,u,~,R_square] = ...
            miaovonmisesfit_ori(deg2rad(orientations),tuning(iCell,:));
        y_fit(:,iCell) = b+R.*exp(k.*(cos(2.*(deg2rad(orisSmooth)-u))-1));
    end    
    vonMisesFits_delays{i} = y_fit;
end
%%
nOrisSmooth = length(orisSmooth);

[~,oriPref_fit] = max(vonMisesFit,[],1);
centeredFits = nan(nOrisSmooth,nCells);
for i = 1:nCells
    nShift = (oriPref_fit(i) - tuningCenterFit)*-1;
    centeredFits(:,i) = circshift(vonMisesFit(:,i),nShift);
end
    
centeredFits_delays = cell(1,nDelay);
centeredResp_delays = cell(1,nDelay);
for i = 1:nDelay
    tuning_fit = vonMisesFits_delays{i};
    [~,oriPref_fit] = max(tuning_fit,[],1);
    tuning_resp = visDelayTuning{i};
    [~,oriPref_resp] = max(tuning_resp,[],2);
    cf = nan(nOrisSmooth,nCells);
    cr = nan(nOri,nCells);
    for icell = 1:nCells
        nShift = (oriPref_fit(icell) - tuningCenterFit)*-1;
        cf(:,icell) = circshift(tuning_fit(:,icell),nShift);
                
        nShift = (oriPref_resp(icell) - tuningCenterResp)*-1;
        cr(:,icell) = circshift(tuning_resp(icell,:),nShift);
    end
    centeredFits_delays{i} = cf;
    centeredResp_delays{i} = cr;
end

avgFits = cell2mat(...
    cellfun(@(x) mean(x(:,isTuned),2),centeredFits_delays,'unif',0));
avgFitsErr = cell2mat(...
    cellfun(@(x) ste(x(:,isTuned),2),centeredFits_delays,'unif',0));
avgResp = cell2mat(...
    cellfun(@(x) mean(x(:,isTuned),2),centeredResp_delays,'unif',0));
avgRespErr = cell2mat(...
    cellfun(@(x) ste(x(:,isTuned),2),centeredResp_delays,'unif',0)); 
avgTrials = round(mean(nTrialsEaDelay(:)));

figure
suptitle({sprintf('N Trials Per Condition (avg) = %s',num2str(avgTrials));...
    sprintf('Tuned Cells (%s)',num2str(sum(isTuned)))})
subplot 121
plot(orisSmooth,avgFits)
figXAxis([],'Orientation (deg)',[-10 190],orientations,orientations)
figYAxis([],'dF/F',[-0.005 0.3])
figAxForm
title('Fits')
subplot 122
x = orientations;
for i  = 1:nDelay
    y = avgResp(:,i);
    yerr = avgRespErr(:,i);
    errorbar(x,y,yerr)
    hold on
end
legend(strsplit(num2str(delays)))
figXAxis([],'Orientation (deg)',[-10 190],orientations,orientations)
figYAxis([],'dF/F',[-0.005 0.3])
figAxForm
title('Responses')

%%
cellInd = isTuned & isResponsive';

[visOnlyFitPeak, visOnlyPref] = max(vonMisesFit,[],1);
[~, pairedPref] = max(vonMisesFits_delays{1},[],1);
pairedFitPeak = nan(1,nCells);
for i = 1:nCells
    pairedFitPeak(i) = vonMisesFits_delays{1}(visOnlyPref(i),i);
end

respLim = [-0.005 0.6];
figure
suptitle('Tuned Responsive Cells')
subplot 121
x = visOnlyFitPeak(cellInd);
y = pairedFitPeak(cellInd);
scatter(x,y)
hold on
errorbar(mean(x),mean(y),ste(y,2),ste(y,2),ste(x,2),ste(x,2),'ro')
plot(respLim,respLim,'k--')
figXAxis([],'Visual Only',respLim)
figYAxis([],'Visual + Auditory',respLim)
figAxForm
title('Pref Ori Resp (ori from vis only condition)')
subplot 122
scatter(visOnlyPref(cellInd)-1,pairedPref(cellInd)-1)
hold on
plot(0:180,0:180,'k--')
figXAxis([],'Visual Only',[0 180])
figYAxis([],'Visual + Auditory',[0 180])
figAxForm
title('Fit Peak Orientation')

fitPeak_delays = nan(nDelay,nCells);
for i = 1:nCells
    fitPeak_delays(:,i) = cellfun(@(x) x(visOnlyPref(i),i),...
        vonMisesFits_delays);
end
avgPeakResp_delays = mean(fitPeak_delays(:,cellInd),2);
errPeakResp_delays = ste(fitPeak_delays(:,cellInd),2);

figure
bar(avgPeakResp_delays)
hold on
errorbar(1:nDelay,avgPeakResp_delays,errPeakResp_delays,'k.')
figXAxis([],'Aud-Vis Delay (ms)',[0 nDelay+1],1:nDelay,delays)
figYAxis([],'dF/F',[])
figAxForm([],0)
title({'Response to vis only pref ori';...
    sprintf('Responsive & Tuned (%s)',num2str(sum(cellInd)))})
%% time-course

[~,oriPref_delay] = cellfun(@(x) max(x,[],2),visDelayTuning,'unif',0);

visPrefTC_delay = cell(1,nDelay);
for i = 1:nDelay
        tc = visTCTuning{i};
        op = oriPref_delay{i};
        ptc = nan(size(tc,1),nCells);
    for icell = 1:nCells
        ptc(:,icell) = tc(:,icell,op(icell));
    end
    visPrefTC_delay{i} = ptc;
end

tcLengthFr = nBaselineFr*4;
avgTC_delay = cellfun(@(x) mean(x(1:tcLengthFr,isTuned),2),...
    visPrefTC_delay,'unif',0);
avgTCErr_delay = cellfun(@(x) ste(x(1:tcLengthFr,isTuned),2),...
    visPrefTC_delay,'unif',0);

colors = brewermap(nDelay+2,'Blues');

tt = round(((1:tcLengthFr) - (nBaselineFr+params.nFramesVisDelay))...
    ./params.frameRate,2,'significant');

figure;
leg = [];
for i = 1:nDelay
    h = shadedErrorBar(tt,avgTC_delay{i},avgTCErr_delay{i},'k');
    h.mainLine.Color = colors(i+2,:);
    h.mainLine.LineWidth = 1;
    h.edge(1).Color = [1 1 1];
    h.edge(2).Color = [1 1 1];
    h.patch.FaceColor = [0.9 0.9 0.9];
    hold on
    leg(i) = h.mainLine;
end
legend(leg,strsplit(num2str(delays)))