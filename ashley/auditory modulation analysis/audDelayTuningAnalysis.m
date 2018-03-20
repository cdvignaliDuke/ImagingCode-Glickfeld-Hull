ds = 'audDelay_V1_EMX';
eval(ds)
rc = behavConstsAV;

fn = fullfile(rc.ashleyAnalysis,'Expt Summaries',ds);
load(fullfile(fn,'dataSummary'))

doLoadPreviousTuning = true;

basewinS = 0.1;
respwinS = 0.1;

minVisOnlyMs = 8000;
nBoot = 1000;

orisSmooth = 0:180;
tuningCenterFit = 90;
fitReliabilityCutoff =85;

nBaselineFr = params.nBaselineMs.*params.frameRate./1000;
basewin = nBaselineFr-round(basewinS*params.frameRate):nBaselineFr;
respwin = (nBaselineFr+params.nFramesVisDelay):...
    (nBaselineFr+params.nFramesVisDelay+(respwinS.*params.frameRate));


im = 1;
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
for i = 1:nOri
    ind = tOris == orientations(i) & tDelay >= 6000;
    visTuning(:,i) = mean(visResp(:,ind),2);
    visTuningErr(:,i) = ste(visResp(:,ind),2);
    nTrialsPerOri(i) = sum(ind);
end

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
isTuned = fitReliability > fitReliabilityCutoff;

delays = unique(tDelay);
nDelay = length(delays);

visDelayTuning = cell(1,nDelay);
visDelayTuningErr = cell(1,nDelay);
visDelayTuningTTest = cell(1,nDelay);
visDelayRespTTest = cell(1,nDelay);
visTCTuning = cell(1,nDelay); 
visTCTuningErr = cell(1,nDelay);
nTrialsEaDelay = zeros(nOri,nDelay);
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
isResponsive = cellfun(@(x,y) x | y,...
                cellfun(@(x) sum(x,2) > 0,visDelayTuningTTest,'unif',0),...
                cellfun(@logical,visDelayRespTTest,'unif',0),'unif',0);

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
figYAxis([],'dF/F',[-0.005 0.1])
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
figYAxis([],'dF/F',[-0.005 0.1])
figAxForm
title('Responses')

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