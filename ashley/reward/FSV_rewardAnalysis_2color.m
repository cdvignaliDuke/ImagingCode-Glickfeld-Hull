clear all
close all
ds = 'FSAV_V1_SOM_temp';
cellsOrDendrites = 1;
doLoadPreviousAnalysis = false;
%%
rc = behavConstsAV;
imgParams_FSAV
bxParams_FSAV_attnV1ms

eval(ds)
titleStr = ds(6:end);
mice = unique({expt.SubNum});
mouse_str = ['i' strjoin(mice,'_i')];

load(fullfile(rc.caOutputDir,ds,...
    [mouse_str '_rewardStruct_cells' ds(5:end) '.mat']));
if strcmp(ds,'FSAV_attentionV1') | strcmp(ds,'FSAV_V1_SOM_temp')
    fnout = fullfile(rc.ashleyAnalysis, 'Expt summaries','reward','behavior',...
        [titleStr '_']); 
    bxExpt = true;
elseif strcmp(ds,'FSAV_V1_SOM_naive_temp')
    fnout = fullfile(rc.ashleyAnalysis, 'Expt summaries','adaptation','naive',...
        [titleStr '_']);     
    bxExpt = false;
end


%%
alignStim = 1;
alignReward = 2;

rewwin = 33:47;

nav = 2;
nBaselineFr = mouse(1).expt(1).info.preAlignFrames;
nFrames1s = frameRateHz;
nexp = size(expt,2);
nCycles = 7;
% lateCycles = 5:nCycles;
lateWinFr = (45:88)+nBaselineFr;
firstWinFr = (3:44)+nBaselineFr;
minTargetRT = (nVisDelayFr_target+respwin_target(1)-nBaselineFr)./frameRateHz*1000;

oriBinSize = 45;
orientations = 0:oriBinSize:(180-oriBinSize);
oriBinEdges = [0, (oriBinSize/2):oriBinSize:(180-(oriBinSize/2)), 180];
nOri = length(orientations);

nMovWin = 15;
movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelFr_target = 30:(30+nMovWin-1);
% movWinLabelMs = 

% minCellN_SIFRmatch = 36;

trOutType = {'h';'m';'fa';'cr'};
trOutTypeName = {'H-All';'H-HT';'H-ET';'M-All';'M-HT';'M-ET';'FA';'CR'};
%% pool experiment data
load([fnout 'adaptAnalysis'],'respCellsExpt')
allRespCellsExpt = struct;
for iexp = 1:size(respCellsExpt,2)
    for itag = 1:2
        allRespCellsExpt(iexp).tag(itag).ind = ...
            respCellsExpt(iexp).tag(itag).firstRespCells | ...
            respCellsExpt(iexp).tag(itag).lateCycRespCells | ...
            respCellsExpt(iexp).tag(itag).targetRespCells ;
    end
end
    rewDataExpt = struct;
%     oriTuningExpt = struct;
    for imouse = 1:size(mouse,2)
        for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end

            d = mouse(imouse).expt(iexp);
            exptID = strcmp({expt.SubNum},mouse(imouse).mouse_name) & strcmp({expt.date},d.date);

            cycLengthFr = d.info.cycTimeFrames;
            nCycLongTC = ceil(longTrialLengthFr./cycLengthFr);

            maxCycles = max(cat(2,d.tag(1).av(visualTrials).align(alignReward).nCycles,...
                d.tag(1).av(auditoryTrials).align(alignReward).nCycles));

            rewDataExpt(exptN).exptName = [mouse(imouse).mouse_name '-' d.date];
            rewDataExpt(exptN).exptCycLengthFr = cycLengthFr;
           itag = 1;
                dd = d.tag(itag).av(visualTrials).align(alignReward);
                cind = logical(allRespCellsExpt(iexp).tag(itag).ind);
                hits = strcmp(dd.outcome,'success');
                misses = strcmp(dd.outcome,'ignore');
                
%                 [nfr,nc,~] = size(dd.respTC);
                orientations = unique(dd.ori);
                nori = length(orientations);
                rewTC = cell(2,nori);
                for iori = 1:nori
                    ind = dd.ori == orientations(iori) & hits;
                    if sum(ind) >= minTrN
                        rewTC{1,iori} = dd.respTC(:,cind,ind);
                    end
                    ind = dd.ori == orientations(iori) & misses;
                    if sum(ind) >= minTrN
                        rewTC{2,iori} = dd.respTC(:,cind,ind);
                    end
                end
                
                hr = squeeze(mean(dd.respTC(rewwin,cind,hits & dd.ori < 30),1));
                mr = squeeze(mean(dd.respTC(rewwin,cind,misses & dd.ori < 30),1));
                rt = mean(dd.reactTime(hits & dd.ori < 30))./1000;
                misswin = (0:14)+round(((nBaselineFr+3+(rt*frameRateHz))/2));
                dm = d.tag(itag).av(visualTrials).align(alignStim);
                mr_stimalign = squeeze(mean(dm.respTC(misswin,cind,...
                    strcmp(dm.outcome,'ignore') & dm.ori <30)));
                auroc_hm1 = nan(1,size(hr,1));
                auroc_hm2 = nan(1,size(hr,1));
                for icell = 1:size(hr,1)
                    auroc_hm1(icell) = roc_gh(mr_stimalign(icell,:),hr(icell,:));
                    auroc_hm2(icell) = roc_gh(mr(icell,:),hr(icell,:));
                end
                
                rewDataExpt(exptN).tag(itag).name = d.tag(itag).name;
                rewDataExpt(exptN).tag(itag).orientations = orientations;
                rewDataExpt(exptN).tag(itag).tc = rewTC;
                rewDataExpt(exptN).tag(itag).aurocReleaseAlign = auroc_hm2;
                rewDataExpt(exptN).tag(itag).aurocStimAlign = auroc_hm1;
           
        end
    end

rewResp = struct;
rewResp.tc = cell(1,2);
rewResp.aurocReleaseAlign = [];
rewResp.aurocStimAlign = [];
for iexp = 1:nexp
    ind = find(rewDataExpt(iexp).tag(1).orientations < 30 &...
        cellfun(@(x) ~isempty(x),rewDataExpt(iexp).tag(1).tc(1,:)));
    tc_hit = [];
    for i = 1:length(ind)
        tc_hit = cat(3,tc_hit,rewDataExpt(iexp).tag(1).tc{1,ind(i)});
    end
    if iexp == 1
        rewResp.tc{1,1} = mean(tc_hit,3);
    else
        rewResp.tc{1,1} = cat(2,rewResp.tc{1,1},mean(tc_hit,3));
    end
    ind = find(rewDataExpt(iexp).tag(1).orientations < 30 &...
        cellfun(@(x) ~isempty(x),rewDataExpt(iexp).tag(1).tc(2,:)));
    tc_miss = [];
    for i = 1:length(ind)
        tc_miss = cat(3,tc_miss,rewDataExpt(iexp).tag(1).tc{2,ind(i)});
    end
    if iexp == 1
        rewResp.tc{1,2} = mean(tc_miss,3);
    else
        rewResp.tc{1,2} = cat(2,rewResp.tc{1,2},mean(tc_miss,3));
    end
    
    rewResp.aurocReleaseAlign = cat(2,rewResp.aurocReleaseAlign,...
        rewDataExpt(iexp).tag(1).aurocReleaseAlign);
    rewResp.aurocStimAlign = cat(2,rewResp.aurocStimAlign,...
        rewDataExpt(iexp).tag(1).aurocStimAlign);    
end
    
%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.01 0.12];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];
siOriLim = [-3 3];
oriBarLim_win = [0 0.08];
oriBarLim_resp = [0 0.04];
oriLim_taskResp = [-0.005 0.035];
oriNLim = [0 120];
oriTCLim = [-0.005 0.08];
targetTCLim = [-0.015 0.08];
outTCLim = [-0.005 0.04];
firstTCLim = [-0.005 0.04];
adaptLim = [0 1];
suppTCLim = [-0.05 0.005];
suppScatLim_win = [-0.2 0.1];
suppScatLim_cyc = [-0.015 0.015];

nTCFrames = size(rewResp.tc{1},1);
tcStartFrame = 26;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = size(antiAnalysis.tag(1).longTC,1);
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

ttMs_rew = ((1:nTCFrames)-nBaselineFr-3).*frameRateHz;
ttFrInd_rew = ismember(round(ttMs_rew,3,'significant'),-1000:500:2000);

nFr_cyc = size(antiAnalysis.tag(itag).cycTC{1,1},1);
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
% tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

lateWinTT = ([lateWinFr(1) lateWinFr(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT = ([respwin(1) respwin(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);
respWinTT_target = (...
    [respwin_target(1) respwin_target(end)] - (nBaselineFr+nVisDelayFr_target))...
    .*(1000/frameRateHz);
baseWinTT = (...
    [basewin_0(1) basewin_0(end)] - (nBaselineFr+nVisDelayFr))...
    .*(1000/frameRateHz);

movWinLabelFr = 30:(30+nMovWin-1);
movWinLabelMs = (movWinLabelFr - (nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);

weightLim = [-3 4.2];
binnedWeightLim = [-0.4 0.4];
weightLimSum = [-0.8 0.8];
siLimSum = [-0.5 2.5];

% lateCycRespAll = mean(antiAnalysis.lateCycTC{1}(respwin,:),1);
%%

hmLim = [-0.3 0.3];
figure
suptitle('Release Aligned: All Resp. Cells, Hard Targets')
colormap(brewermap([],'*RdBu'));
itag = 1;
subplot(2,3,1)
resp = mean(rewResp.tc{1}((nBaselineFr+1):end,:),1);
[~,respSortInd] = sort(resp);
hm = flipud(rewResp.tc{1}(:,respSortInd)');
imagesc(hm)
hold on
figXAxis([],'Time from Release (ms)',[],...
    find(ttFrInd_rew),ttMs_rew(ttFrInd_rew))
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title('Hits - dF/F')

subplot(2,3,2)
% resp = mean(rewResp.tc{2}((nBaselineFr+1):end,:),1);
% [~,respSortInd] = sort(resp);
hm = flipud(rewResp.tc{2}(:,respSortInd)');
imagesc(hm)
hold on
figXAxis([],'Time from Release (ms)',[],...
    find(ttFrInd_rew),ttMs_rew(ttFrInd_rew))
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title('Misses - dF/F')

subplot(2,3,4)
histogram(rewResp.aurocReleaseAlign,0:0.1:1)
figXAxis([],'auROC',[0 1])
figYAxis([],'N Cells',[0 25])
figAxForm
title('Hit vs. Release Aligned Miss')

subplot(2,3,5)
histogram(rewResp.aurocStimAlign,0:0.1:1)
figXAxis([],'auROC',[0 1])
figYAxis([],'N Cells',[0 25])
figAxForm
title('Hit vs. Stim Aligned Miss')

subplot(2,3,6)
plot(rewResp.aurocReleaseAlign,rewResp.aurocStimAlign,'.','MarkerSize',10)
hold on
plot(0:1,0:1,'k--')
figXAxis([],'Release',[0 1])
figYAxis([],'Stim',[])
figAxForm
title('Hit vs. Miss auROC')

print([fnout 'rewardAuROC_HM'],'-dpdf','-fillpage')
%%
    if bxExpt
        subplot(2,2,itag+2)
        lateWinTC = antiAnalysis_passive.tag(itag).longTC(:,ind);
        hm = flipud(lateWinTC(:,respSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        hold on
        if strcmp(ds,'FSAV_attentionV1')
            exCellInd = [exampleCell_1,exampleCell_2];
            exCellMat = zeros(1,length(cellInfo.firstRespCells));
            exCellMat(exCellInd) = 1;
            exCellSortInd = find(flip(exCellMat(respSortInd)));
            hline(exCellSortInd,'k-')
        end
        figXAxis([],'Time from Start (fr)',[],find(ttFrInd_rew),ttMs_rew(ttFrInd_rew))
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('Passive, %s Cells',antiAnalysis.tag(itag).name))
    end
end
print([fnout 'heatmapsAllCells'],'-dpdf','-fillpage')

figure
suptitle('Easy Targets (>60 deg): All Cells')
colormap(brewermap([],'*RdBu'));
for itag = 1:2
    subplot(2,2,itag)
    ind = logical(cellInfo.tag(itag).targetRespCells);
%     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
%     ind = cellInfo.tag(itag).lateWinRespCells | cellInfo.tag(itag).lateWinSuppCells;
    tc = tarAnalysis.tag(itag).easyTC(:,ind);
    resp = mean(tc(respwin,:),1);
    [~,respSortInd] = sort(resp);
    hm = flipud(tc(:,respSortInd)');
    imagesc(hm(:,tcStartFrame:end))
    hold on
    figXAxis([],'Time from Start (fr)',[],ttLabelFr_target,ttLabelFr_target)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title(sprintf('Behaving, %s Cells',antiAnalysis.tag(itag).name))
    
    if bxExpt
        subplot(2,2,itag+2)
        tc = tarAnalysis_passive.tag(itag).easyTC(:,ind);
        hm = flipud(tc(:,respSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        hold on
        figXAxis([],'Time from Start (fr)',[],ttLabelFr_target,ttLabelFr_target)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('Passive, %s Cells',antiAnalysis.tag(itag).name))
    end
end

hmLim = [-0.1 0.1];
figure
suptitle('Late Resp. Cells')
colormap(brewermap([],'*RdBu'));
for itag = 1:2
    subplot(2,2,itag)
    ind = logical(cellInfo.tag(itag).lateCycRespCells);
%     ind = logical(cellInfo.tag(itag).lateWinSuppCells);
%     ind = cellInfo.tag(itag).lateWinRespCells | cellInfo.tag(itag).lateWinSuppCells;
    tc = antiAnalysis.tag(itag).lateCycTC(:,ind);
    resp = mean(tc(respwin,:),1);
    [~,respSortInd] = sort(resp);
    hm = flipud(tc(:,respSortInd)');
    imagesc(hm(:,tcStartFrame:end))
    hold on
    figXAxis([],'Time from Start (fr)',[],ttLabelFr_target,ttLabelFr_target)
    figYAxis([],'Cell #',[])
    figAxForm
    colorbar
    caxis(hmLim)
    title(sprintf('Behaving, %s Cells',antiAnalysis.tag(itag).name))
    
    if bxExpt
        subplot(2,2,itag+2)
        tc = antiAnalysis_passive.tag(itag).lateCycTC(:,ind);
        hm = flipud(tc(:,respSortInd)');
        imagesc(hm(:,tcStartFrame:end))
        hold on
        figXAxis([],'Time from Start (fr)',[],ttLabelFr_target,ttLabelFr_target)
        figYAxis([],'Cell #',[])
        figAxForm
        colorbar
        caxis(hmLim)
        title(sprintf('Passive, %s Cells',antiAnalysis.tag(itag).name))
    end
end

%% cluster analysis
allCells = [];
allCells_pass = [];
cellTag = [];
for itag = 1:2
    ind = cellInfo.tag(itag).lateWinRespCells | ...
        cellInfo.tag(itag).lateWinSuppCells | ...
        cellInfo.tag(itag).firstRespCells | ...
        cellInfo.tag(itag).lateCycRespCells;
    lateWinTC = antiAnalysis.tag(itag).longTC(:,ind);
    lateWinResp = mean(lateWinTC(lateWinFr,:),1);
    [~,respSortInd] = sort(lateWinResp);
    lateWinTC_pass = antiAnalysis_passive.tag(itag).longTC(:,ind);
    allCells = cat(2,allCells,lateWinTC(:,respSortInd));
    allCells_pass = cat(2,allCells_pass,lateWinTC_pass(:,respSortInd));
    cellTag = cat(2,cellTag,repmat({cellInfo.tag(itag).name},1,sum(ind)));
end

nc = size(allCells,2);
allCells_norm = allCells./max(allCells,[],1);
allCells_pass_norm = allCells_pass./max(allCells,[],1);

% allCells_pca = pca(allCells_norm);

nCluster = 3;
[clustID, clustCentroid] = kmeans(allCells_norm',nCluster,'MaxIter',1000000,'Replicates',1000);
[~,clustID_pass] = pdist2(clustCentroid,allCells_pass_norm','euclidean','smallest',1);

figure
subplot 221
histogram(clustID,0:(nCluster+1),'Normalization','probability')
figXAxis([],'Cluster #',[0 (nCluster+2)])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
subplot 222
histogram(clustID_pass,0:(nCluster+1),'Normalization','probability')
figXAxis([],'Cluster #',[0 (nCluster+2)])
figYAxis([],'Fraction of Cells',[0 1])
figAxForm
subplot 223
for ic = 1:nCluster
    cInd = clustID == ic;
    hold on
    plot(ones(1,sum(cInd)).*ic,find(cInd),'.','MarkerSize',10);
end
tagInd = strcmp(cellTag,'SOM+');
plot(ones(1,sum(tagInd)).*(nCluster+1),find(tagInd),'k.','MarkerSize',10)
figXAxis([],'Cluster #',[0 (nCluster+2)],1:(nCluster+1),...
    cat(2,cellfun(@num2str,num2cell(1:nCluster),'unif',0),{'SOM+'}))
figYAxis([],'Fraction of Cells',[0 nc])
figAxForm
subplot 224
for ic = 1:nCluster
    cInd = clustID_pass == ic;
    hold on
    plot(ones(1,sum(cInd)).*ic,find(cInd),'.','MarkerSize',10);
end
figXAxis([],'Cluster #',[0 (nCluster+1)])
figYAxis([],'Fraction of Cells',[0 nc])
figAxForm
print([fnout 'clusterAnalysisBxPass_clusterID'],'-dpdf','-fillpage')

clusterTC = cell(1,nCluster);
clusterTC_pass = cell(1,nCluster);
for ic = 1:nCluster
    cInd = clustID == ic;
    clusterTC{ic} = allCells(:,cInd);
    cInd = clustID_pass == ic;
    clusterTC_pass{ic} = allCells_pass(:,cInd);    
end

figure
for ic = 1:nCluster
    subplot(2,nCluster,ic)
    y = mean(clusterTC{ic}(tcStartFrame:end,:),2);
    yerr = ste(clusterTC{ic}(tcStartFrame:end,:),2);
    if all(~isnan(yerr))
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    else
        plot(tt_longTC,y,'k')
    end
    figXAxis([],'Time from start (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long)
    figYAxis([],'dF/F',[-0.1 0.1])
    figAxForm
    title(sprintf('Cluster %s, n=%s',num2str(ic),num2str(size(clusterTC{ic},2))))
    
    subplot(2,nCluster,ic+nCluster)
    y = mean(clusterTC_pass{ic}(tcStartFrame:end,:),2);
    yerr = ste(clusterTC_pass{ic}(tcStartFrame:end,:),2);
    if all(~isnan(yerr))
        shadedErrorBar_chooseColor(tt_longTC,y,yerr,[0 0 0]);
    else
        plot(tt_longTC,y,'k')
    end
    figXAxis([],'Time from start (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long)
    figYAxis([],'dF/F',[-0.1 0.1])
    figAxForm
    title(sprintf('Cluster %s, n=%s',num2str(ic),num2str(size(clusterTC_pass{ic},2))))
end
print([fnout 'clusterAnalysisBxPass_clusterTC'],'-dpdf','-fillpage')

switchClustID = clustID'-clustID_pass ~= 0;

figure
subplot 121
y = sum(switchClustID & ~tagInd)./sum(~tagInd);
plot(1,y,'.','MarkerSize',20)
hold on
y = sum(switchClustID & tagInd)./sum(tagInd);
plot(2,y,'.','MarkerSize',20);
figXAxis([],'',[0 3],1:2,{'SOM-','SOM+'})
figYAxis([],'Fraction Switch Cluster',[0 0.2])
figAxForm
subplot 122
hold on
for ic = 1:nCluster
    cInd = clustID' == ic & switchClustID;
    n = histcounts(clustID_pass(cInd),1:nCluster+1);
    scatter(ic.*(ones(1,nCluster)),1:nCluster,(n*5)+1)
end
figXAxis([],'Behavior Cluster #',[0 nCluster+1])
figYAxis([],'Passive Cluster #',[0 nCluster+1])
figAxForm
title(sprintf('Cluster Switch Cells, n=%s',num2str(sum(switchClustID))))
print([fnout 'clusterAnalysisBxPass_nClusterSwitch'],'-dpdf','-fillpage')

%%
    figure
    suptitle('Target Resp.')
for itag = 1:2
    ind = logical(cellInfo.tag(itag).targetRespCells);
    tarResp_bx = mean(tarAnalysis.tag(itag).easyTC(respwin,:),1) - mean(tarAnalysis.tag(itag).easyTC(basewin_0,:),1);
    tarResp_pass = mean(tarAnalysis_passive.tag(itag).easyTC(respwin,:),1) - mean(tarAnalysis_passive.tag(itag).easyTC(basewin_0,:),1);

    tarRespLim = [-0.01 0.12];
%     for icyc = 1:nCycles
        subplot(1,2,itag)
        x = tarResp_bx(ind);
        y = tarResp_pass(ind);
        plot(x,y,'.','MarkerSize',10)
        hold on
        plot(tarRespLim,tarRespLim,'k--')
        figXAxis([],'Behavior',tarRespLim)
        figYAxis([],'Passive',tarRespLim)
        figAxForm
        title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
%     end
end

figure
    suptitle('Late Resp.')
for itag = 1:2
    ind = logical(cellInfo.tag(itag).lateCycRespCells);
    lateResp_bx = mean(antiAnalysis.tag(itag).lateCycTC(respwin,:),1) - mean(antiAnalysis.tag(itag).lateCycTC(basewin_0,:),1);
    lateResp_pass = mean(antiAnalysis_passive.tag(itag).lateCycTC(respwin,:),1) - mean(antiAnalysis_passive.tag(itag).lateCycTC(basewin_0,:),1);

    lateRespLim = [-0.01 0.1];
%     for icyc = 1:nCycles
        subplot(1,2,itag)
        x = lateResp_bx(ind);
        y = lateResp_pass(ind);
        plot(x,y,'.','MarkerSize',10)
        hold on
        plot(lateRespLim,lateRespLim,'k--')
        figXAxis([],'Behavior',lateRespLim)
        figYAxis([],'Passive',lateRespLim)
        figAxForm
        title(sprintf('%s Cells',antiAnalysis.tag(itag).name))
%     end
end