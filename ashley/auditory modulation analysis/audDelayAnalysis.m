clear all
close all
ds = '_V1';
%%
rc = behavConstsAV;
eval(['awData_audMod' ds])

%% load data
fn = fullfile(rc.ashleyAnalysis,'Expt Summaries',['awData_audMod' ds]);
load(fullfile(fn,['awData_audMod' ds '_CaSummary']),'ms');

%% some params

preFrames = ms(1).pre_event_frames;
ntrtype = length(ms(1).delayTrials);
trType_str = ms(1).audDelayType;
frRateHz = ms(1).frRateHz;
oris = [];
for iexp = 1:size(ms,2)
    msOris = ms(iexp).orientations;
    oris = cat(2,oris,msOris);
end   
oris = unique(oris);
nori = length(oris);
  
pre_win = preFrames-round(.25*frRateHz)+1:preFrames;
resp_win = preFrames + round(.1*frRateHz):preFrames + round(.35*frRateHz);

c = brewermap(ntrtype+1,'*Blues');
delay_colors(1:ntrtype-1,:) = c(2:end-1,:);
delay_colors(ntrtype,:) = [0 0 0];
%% combine datasets
delayTCOri = cell(noris,ntrtype);
delayTCOri_err = cell(noris,ntrtype);
delayOri_ttest = cell(noris,ntrtype);
for iexp = 1:size(ms,2)
    for itype = 1:ntrtype
        dt = ms(iexp).delayTrials{itype};
        for iori = 1:nori
            ind = ms(iexp).trOri{itype} == oris(iori);
            delayTCOri{iori,itype} = mean(dt(:,:,ind),3);
            delayTCOri_err{iori,itype} = ste(dt(:,:,ind),3);
            pt = mean(dt(pre_win,:,ind),1);
            rt = mean(dt(resp_win,:,ind),1);
            delayOri_ttest{iori,itype} = ttest(rt,pt,'dim',3,'tail','right');
        end
    end
end

delayRespOri = cellfun(@(x) mean(x(resp_win,:),1) - mean(x(pre_win,:),1), delayTCOri,'unif',0);
delayRespOri_err = cellfun(@(x) mean(x(resp_win,:),1) - mean(x(pre_win,:),1), delayTCOri_err,'unif',0);
%% cell orientation tuning
nc = size(delayTCOri{1},2);
delayTuning = cell(1,ntrtype);
delayTuning_err = cell(1,ntrtype);
delayTuningTtest = cell(1,ntrtype);
for itype = 1:ntrtype
    delayTuning{itype} = cell2mat(delayRespOri(:,itype));
    delayTuning_err{itype} = cell2mat(delayRespOri(:,itype));
    delayTuningTtest{itype} = cell2mat(delayOri_ttest(:,itype));
end
[~, delayOriMax_ind] = cellfun(@(x) max(x,[],1),delayTuning,'unif',0);

% for vis only data, find if max response is also significant response
visOnlyMaxInd = delayOriMax_ind{strcmp(trType_str,'vis only')};
visOnlyTtest = delayTuningTtest{strcmp(trType_str,'vis only')};
visOnlyMaxSN = zeros(1,nc);
for icell = 1:nc
    tt = visOnlyTtest(:,icell);
   visOnlyMaxSN(icell) = tt(visOnlyMaxInd(icell));
end

visOnlyMaxInd_sub1 = visOnlyMaxInd-1;
visOnlyMaxInd_sub1(visOnlyMaxInd_sub1 == 0) = nori;
%% plot example cells orientation curves for each delay
exCells = randsample(find(visOnlyMaxSN),16);
exC = exCells(1);

% for icell = 1:16
%     exC = exCells(icell);
figure;
for itype = 1:ntrtype
    tc = delayTuning{itype}(:,exC);
    err = delayTuning_err{itype}(:,exC);
    hold on
    h = errorbar(oris,tc,err,'ko-');
    h.Color = delay_colors(itype,:);
    h.LineWidth = 1;
    h.MarkerFaceColor = delay_colors(itype,:);
end
legend(trType_str,'location','northeastoutside')
figXAxis([],'orientation (deg)',[oris(1)-10 oris(end)+10],oris,oris)
figYAxis([],'dF/F',[])
figAxForm([])
hold on
plot(oris(visOnlyMaxInd(exC)),tc(visOnlyMaxInd(exC))+0.1,'k*')
title(['cell #' num2str(exC)])
% end

print(fullfile(fn,'exCell_tuning'),'-dpdf','-fillpage');

%% plot example cell time-course for each delay, multiple orientations
trLenFr = size(delayTCOri{1,1},1);
trMs = round((-preFrames+1:trLenFr-preFrames)/frRateHz*1000);
exCells = randsample(find(visOnlyMaxSN),16);
exC = exCells(1);
figure;
suptitle(['example cell #' num2str(exC)])
for iori = 1:nori
    tc1 = delayTCOri{iori,1}(:,exC);
    err1 = delayTCOri_err{iori,1}(:,exC);
    tc2 = delayTCOri{iori,5}(:,exC);
    err2 = delayTCOri_err{iori,5}(:,exC);
    subplot(1,nori,iori)
    h = shadedErrorBar(trMs,tc1,err1,'m-');
    leg(1) = h.mainLine;
    hold on
    h = shadedErrorBar(trMs,tc2,err2,'k-');
    leg(2) = h.mainLine;
    hold on
    vline(0,'k--')
    figXAxis([],'time from vis stim on (ms)',[-250 1000])
    figYAxis([],'dF/F',[-0.05 0.15])
    figAxForm([])
    title([num2str(oris(iori)) ' deg'])
    if iori == 1
        legend(leg,{'vis+aud';'vis only'},'location','northwest');
    end
end

print(fullfile(fn,'exCell_tc_tuning'),'-dpdf','-fillpage');
%% plot scatter of responses to each delay type compared to vis only 
% (prefferred orienation)

% get response in pref ori, only for signif cells
visOnlyResp_pref = NaN(1,nc);
audOnlyResp_pref = NaN(1,nc);
delayResp_pref = cell(1,4);
for icell = 1:nc
    ind = visOnlyMaxInd(icell);
    if visOnlyMaxSN(icell)
        visOnlyResp_pref(icell) = delayTuning{5}(ind,icell);
        audOnlyResp_pref(icell) = delayTuning{6}(ind,icell);
        for idelay = 1:4
            delayResp_pref{idelay}(icell) = delayTuning{idelay}(ind,icell);
        end
    else
        visOnlyResp_pref(icell) = NaN;
        audOnlyResp_pref(icell) = NaN;
        for idelay = 1:4
            delayResp_pref{idelay}(icell) = NaN;
        end
    end
end

%plot scatter
resp_lim = [-0.05 0.15];
figure;
suptitle({'response to preferred orientation'; 'preference is max response on vis only trials, if significant'})
for idelay = 1:4
    subplot(2,2,idelay)
    scatter(visOnlyResp_pref,delayResp_pref{idelay},50,'k.');
    hold on 
    plot(resp_lim,resp_lim,'k--')
    figXAxis([],trType_str{5},resp_lim)
    figYAxis([],[trType_str{idelay} 'ms delay'],resp_lim)
    figAxForm([])
end
print(fullfile(fn,'scatter_maxResp_eaDelay'),'-dpdf','-fillpage');
%% plot scatter of responses to each delay type compared to vis only 
% (prefferred - 1 orienation)

% get response in pref ori, only for signif cells
visOnlyResp_prefsub1 = NaN(1,nc);
audOnlyResp_prefsub1 = NaN(1,nc);
delayResp_prefsub1 = cell(1,4);
for icell = 1:nc
    ind = visOnlyMaxInd_sub1(icell);
    if visOnlyMaxSN(icell)
        visOnlyResp_prefsub1(icell) = delayTuning{5}(ind,icell);
        audOnlyResp_prefsub1(icell) = delayTuning{6}(ind,icell);
        for idelay = 1:4
            delayResp_prefsub1{idelay}(icell) = delayTuning{idelay}(ind,icell);
        end
    else
        visOnlyResp_prefsub1(icell) = NaN;
        audOnlyResp_prefsub1(icell) = NaN;
        for idelay = 1:4
            delayResp_prefsub1{idelay}(icell) = NaN;
        end
    end
end

%plot scatter
resp_lim = [-0.05 0.15];
figure;
suptitle({'response to pref-adjacent orientation'; 'preference is max response on vis only trials, if significant'})
for idelay = 1:4
    subplot(2,2,idelay)
    scatter(visOnlyResp_prefsub1,delayResp_prefsub1{idelay},50,'k.');
    hold on 
    plot(resp_lim,resp_lim,'k--')
    figXAxis([],trType_str{5},resp_lim)
    figYAxis([],[trType_str{idelay} 'ms delay'],resp_lim)
    figAxForm([])
end
print(fullfile(fn,'scatter_maxRespSub1_eaDelay'),'-dpdf','-fillpage');
%% compare prefferred stim and prefsub1 stim
allDelayResp_pref = cell2mat(cat(2,delayResp_pref,visOnlyResp_pref,audOnlyResp_pref)');
allDelayPref = mean(allDelayResp_pref(:,logical(visOnlyMaxSN)),2);
allDelayPref_err = ste(allDelayResp_pref(:,logical(visOnlyMaxSN)),2);

allDelayResp_sub1 = cell2mat(cat(2,delayResp_prefsub1,visOnlyResp_prefsub1,audOnlyResp_prefsub1)');
allDelaySub1 = mean(allDelayResp_sub1(:,logical(visOnlyMaxSN)),2);
allDelaySub1_err = ste(allDelayResp_sub1(:,logical(visOnlyMaxSN)),2);

figure;
h = errorbar(1:ntrtype,allDelayPref,allDelayPref_err,'ko');
h.MarkerFaceColor = 'k';
hold on
h = errorbar(1:ntrtype,allDelaySub1,allDelaySub1_err,'ko');
h.MarkerFaceColor = [0.5 0.5 0.5];
h.Color = [0.5 0.5 0.5];
figXAxis([],'delay from auditory stim',[0 ntrtype+1],1:ntrtype,trType_str)
figYAxis([],'dF/F',[]);
figAxForm([])
legend({'resp 2 pref ori';'resp 2 pref ori-1'},'location','northeastoutside')
title({'response at preferred orienation (responsive cells only)';'error is ste across cells'})

print(fullfile(fn,'avgResp_eaDelay'),'-dpdf','-fillpage');
%% other things to do
% get a good measure of orienation selectivity, i.e. gOSI;
% compare responses to true preferred orientation
