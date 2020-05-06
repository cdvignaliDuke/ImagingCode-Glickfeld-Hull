
nMice = length(mice);
channelColor = {'green','red'};
avName = {'visual','auditory'};

%%

%% plotting params
respTCLim = [-0.005 0.05];
cycTCLim = [-0.005 0.015];
eaCycTCLim = [-0.005 0.05];
cycTCLim_minRespCells = [-0.005 0.025];
scatLim_win = [-0.2 0.6];
scatLim_cyc = [-0.035 0.085];
hmLim = [-0.1 0.1];
exCellTCLim = [-0.02 0.15];
oriRespLim = [-0.05 0.15];
siLim = [-10 10];
siOriLim = [-4 4];
oriBarLim_win = [0 0.08];
oriBarLim_resp = [0 0.04];
oriLim_taskResp = [-0.005 0.035];
oriNLim = [0 120];
oriTCLim = [-0.005 0.08];
targetTCLim = [-0.015 0.08];
outTCLim = [-0.005 0.04];
firstTCLim = [-0.005 0.04];
% adaptLim = [0 1];
suppTCLim = [-0.05 0.005];
suppScatLim_win = [-0.2 0.1];
suppScatLim_cyc = [-0.015 0.015];
cellRespTCLim = [-0.05 0.15];
exCellRespTCLim = [-0.1 0.1];
siBinLim = [-4 4];
stimRespLim = [-0.01 0.05];
avName = {'Vis','Aud'};

tcStartFrame = 26;
tcEndFrame = 45;
cycTCEndTimeMs = 350;
cycTCEndFr = 45;
ttLabel_long = 0:500:2500;
ttLabel_cyc = -200:100:cycTCEndTimeMs;
ttLabel_target = -1000:250:900;
preTargetStimLabel = -700:350:0;
nFr_long = 100;
tt_longTC = ((tcStartFrame:nFr_long)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
ttLabelFr_long = ((ttLabel_long./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_cyc = ((ttLabel_cyc./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr)-tcStartFrame+1);
ttLabelFr_target = ((ttLabel_target./1000)*frameRateHz)+...
    ((nBaselineFr+nVisDelayFr_target)+1);

nFr_cyc = 45;
tt_cycTC = ((tcStartFrame:nFr_cyc)-(nBaselineFr+nVisDelayFr)).*(1000/frameRateHz);
tt_targetTC = ((1:nFr_cyc)-(nBaselineFr+nVisDelayFr_target)).*(1000/frameRateHz);

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

%% grab SOM data
nexp = length(som_expt);
expt_adapt = struct;
expt_sort = nan(1,nexp);
for im = 1:nMice
    ne = size(mouse(im).expt,2);
    if any(strcmp({som_expt.SubNum},mouse(im).mouse_name))
        for iexp = 1:ne
            if isempty(fieldnames(expt_adapt))
                exptN = 1;
            else
                exptN = exptN+1;
            end
            d = mouse(im).expt(iexp);
            expt_adapt(exptN).expt_name = [mouse(im).mouse_name '-' d.date];
            ind = strcmp({expt.SubNum},mouse(im).mouse_name) & ...
                strcmp({expt.date},d.date);
            if expt(ind).indicator{2}(1:4) == 'flex'
                expt_adapt(exptN).indicator = expt(ind).indicator{2}(6:end);
            else
                expt_adapt(exptN).indicator = expt(ind).indicator{2};
            end
            expt_adapt(exptN).isBx = logical(expt(ind).isBehav);
            expt_sort(exptN) = find(ind);

            for itag = 1:2
                if itag == 1
                    tag_name = expt(ind).greenChannelLabel;
                    if strcmp(tag_name,'SOM')
                        somIsGreen = true;
                    else
                        somIsGreen = false;
                    end
                else
                    tag_name = expt(ind).redChannelLabel;
                    if strcmp(tag_name,'SOM')
                        somIsGreen = false;
                    end
                end
                expt_adapt(exptN).tag(itag).name = tag_name;
                if ~isempty(tag_name)
                    for iav = 1:2
                        if somIsGreen                            
                            expt_adapt(exptN).tag(1).av(iav) = d.tag(2).av(iav);
                        else
                            expt_adapt(exptN).tag(itag).av(iav) = d.tag(itag).av(iav);
                        end
                    end
                end

            end
        end
    end
end

indicators = unique({expt_adapt.indicator});

resp_adapt = struct;
for iexp = 1:nexp
    if expt_adapt(iexp).tag(1).name == 'SOM'
        somInd = 1;
    else
        somInd = 2;
    end
    resp_adapt(iexp).isBx = expt_adapt(iexp).isBx;
    resp_adapt(iexp).channel = channelColor{somInd};
    resp_adapt(iexp).indicator = expt_adapt(iexp).indicator;
    if ~isempty(expt_adapt(iexp).tag(somInd).av(1).align(1).respTC)
        for iav = 1:2
            resp_adapt(iexp).av(iav).first = expt_adapt(iexp).tag(somInd).av(iav).align(1).respTC...
                - mean(expt_adapt(iexp).tag(somInd).av(iav).align(1).respTC(basewin,:,:),1);
            resp_adapt(iexp).av(iav).late = expt_adapt(iexp).tag(somInd).av(iav).align(3).respTC...
                - mean(expt_adapt(iexp).tag(somInd).av(iav).align(3).respTC(basewin,:,:),1);
        end
    end
end

cellInfo = struct;
for iexp = 1:nexp
    d_av_first = cat(3,resp_adapt(iexp).av(visualTrials).first,...
        resp_adapt(iexp).av(auditoryTrials).late);
    d_av_late = cat(3,resp_adapt(iexp).av(visualTrials).first,...
        resp_adapt(iexp).av(auditoryTrials).late);
    cellInfo(iexp).responsive_first = ttest(squeeze(mean(d_av_first(respwin,:,:))),...
        squeeze(mean(d_av_first(basewin,:,:))),'dim',2,'tail','right');
    cellInfo(iexp).responsive_first = ttest(squeeze(mean(d_av_late(respwin,:,:))),...
        squeeze(mean(d_av_late(basewin,:,:))),'dim',2,'tail','right');
end
%% plot first to late adaptation across naive an behavior expt
iav = 1;
late_color = [0.4 1 0.4];
early_color = {[0 0 0],[.5 .5 .5]};

setFigParams4Print('portrait')
figure
subplot 321
early_tc = [];
late_tc = [];
ind = [];
for iexp = 1:nexp
    early_tc = cat(2,early_tc,mean(resp_adapt(iexp).av(iav).first,3)-...
        mean(mean(resp_adapt(iexp).av(iav).first(basewin_0,:,:),3)));
    late_tc = cat(2,late_tc,mean(resp_adapt(iexp).av(iav).late,3)-...
        mean(mean(resp_adapt(iexp).av(iav).late(basewin_0,:,:),3)));
    ind = cat(1,ind,cellInfo(iexp).responsive_first);
end
y = mean(early_tc(tcStartFrame:tcEndFrame,ind==1),2);
yerr = ste(early_tc(tcStartFrame:tcEndFrame,ind==1),2);
hold on
shadedErrorBar_chooseColor(tt_cycTC,y,yerr,early_color{1});
y = mean(late_tc(tcStartFrame:tcEndFrame,ind==1),2);
yerr = ste(late_tc(tcStartFrame:tcEndFrame,ind==1),2);
shadedErrorBar_chooseColor(tt_cycTC,y,yerr,late_color);
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',[]) 
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title(sprintf('All Resp. Cells (%s/%s)',...
    num2str(sum(ind)),...
    num2str(length(ind))))

early_tc = cell(1,2);
late_tc = cell(1,2);
ind = cell(1,2);
for iexp = 1:nexp
    if resp_adapt(iexp).isBx == 1
        bxInd = 1;
    else
        bxInd = 2;
    end
    early_tc{bxInd} = cat(2,early_tc{bxInd},mean(resp_adapt(iexp).av(iav).first,3)-...
        mean(mean(resp_adapt(iexp).av(iav).first(basewin_0,:,:),3)));
    late_tc{bxInd} = cat(2,late_tc{bxInd},mean(resp_adapt(iexp).av(iav).late,3)-...
        mean(mean(resp_adapt(iexp).av(iav).late(basewin_0,:,:),3)));
    ind{bxInd} = cat(1,ind{bxInd},cellInfo(iexp).responsive_first);
end
for i = 1:2
    a = mean(late_tc{i}(respwin,ind{i}==1),1);
    a(a<0) = 0;
    b = mean(early_tc{i}(respwin,ind{i}==1),1);
    subplot 323
    y = [mean(b),mean(a)];
    yerr = [ste(b,2),ste(a,2)];
    hold on
    if i == 1
        h = errorbar(1:2,y,yerr,'k.-');
    else
        h = errorbar(1:2,y,yerr,'k.--');
    end
    subplot 325
    hold on
    if i == 1
        h = errorbar(1,mean(a./b),ste(a./b,2),'k.','MarkerSize',20);
    else
        h = errorbar(1,mean(a./b),ste(a./b,2),'ko');
    end
end
subplot 323
figXAxis([],'Time Point',[0 3],1:2,{'Early','Late'})
figYAxis([],'dF/F',[])
figAxForm
legend({'Behav','Naive'},'location','northeast')
subplot 325
figXAxis([],'',[0 2],1,{'All'})
figYAxis([],'Norm. dF/F (Late/Early)',[])
figAxForm

early_tc = cell(1,2);
late_tc = cell(1,2);
ind = cell(1,2);
for iexp = 1:nexp
    if strcmp(resp_adapt(iexp).indicator,'GCaMP6s')
        gcampInd = 1;
    else
        gcampInd = 2;
    end
    early_tc{gcampInd} = cat(2,early_tc{gcampInd},mean(resp_adapt(iexp).av(iav).first,3)-...
        mean(mean(resp_adapt(iexp).av(iav).first(basewin_0,:,:),3)));
    late_tc{gcampInd} = cat(2,late_tc{gcampInd},mean(resp_adapt(iexp).av(iav).late,3)-...
        mean(mean(resp_adapt(iexp).av(iav).late(basewin_0,:,:),3)));
    ind{gcampInd} = cat(1,ind{gcampInd},cellInfo(iexp).responsive_first);
end
for i = 1:2
    subplot 322
    y = mean(early_tc{i}(tcStartFrame:tcEndFrame,ind{i}==1),2);
    yerr = ste(early_tc{i}(tcStartFrame:tcEndFrame,ind{i}==1),2);
    hold on
    shadedErrorBar_chooseColor(tt_cycTC,y,yerr,early_color{i});
    
    subplot 324
    a = mean(late_tc{i}(respwin,ind{i}==1),1);
    a(a<0) = 0;
    b = mean(early_tc{i}(respwin,ind{i}==1),1);
    y = [mean(b),mean(a)];
    yerr = [ste(b,2),ste(a,2)];
    hold on
    h = errorbar(1:2,y,yerr,'.-');
    h.Color = early_color{i};
    
    subplot 326
    hold on    
    h = errorbar(1,mean(a./b),ste(a./b,2),'k.','MarkerSize',20);
    h.Color = early_color{i};
end
subplot 322
figXAxis([],'Time (ms)',[tt_cycTC(1) cycTCEndTimeMs],ttLabel_cyc,ttLabel_cyc)
figYAxis([],'dF/F',[]) 
hline(0,'k:')
vline(respWinTT,'k--')
figAxForm
title('Early Resp')
subplot 324
figXAxis([],'Time Point',[0 3],1:2,indicators)
figYAxis([],'dF/F',[])
figAxForm
legend(indicators,'location','northeast')
subplot 326
figXAxis([],'',[0 2],1,{'All'})
figYAxis([],'Norm. dF/F (Late/Early)',[])
figAxForm


print([fnout 'lateAdaptation'],'-dpdf','-fillpage')

%% each stimulus adapation
cycLengthFr = 11;
stimResp_adapt = struct;
for iexp = 1:nexp
    if expt_adapt(iexp).tag(1).name == 'SOM'
        somInd = 1;
    else
        somInd = 2;
    end
    stimResp_adapt(iexp).isBx = expt_adapt(iexp).isBx;
    stimResp_adapt(iexp).channel = channelColor{somInd};
    stimResp_adapt(iexp).indicator = expt_adapt(iexp).indicator;
    for iav = 1:2
        d = expt_adapt(iexp).tag(somInd).av(iav).align(1);
        ind = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
        cycTC = cell(1,maxCycles);
        for icyc = 1:maxCycles
            tc = d.respTC(:,:,d.nCycles >= icyc & ind);
            cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
            cycTC{1,icyc} = tc(...
                (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
        end
        stimResp_adapt(iexp).av(iav).cycTC = cycTC;
    end
end

%% plot each stimulus adaptation
iav=1;
bxName = {'Behav','Naive'};
hm_lim = [-0.15 0.15];
hm_adapt = cell(1,2);
for iexp = 1:nexp
    if stimResp_adapt(iexp).isBx
        hm_adapt{1} = cat(1,hm_adapt{1},cell2mat(cellfun(@(x)...
            squeeze(mean(mean(x(respwin,cellInfo(iexp).responsive_first==1,:),1)...
            -mean(x(basewin_0,cellInfo(iexp).responsive_first==1,:),1),3))',...
            stimResp_adapt(iexp).av(iav).cycTC,'unif',0)));
    else
        hm_adapt{2} = cat(1,hm_adapt{2},cell2mat(cellfun(@(x)...
            squeeze(mean(mean(x(respwin,cellInfo(iexp).responsive_first==1,:),1)...
            -mean(x(basewin_0,cellInfo(iexp).responsive_first==1,:),1),3))',...
            stimResp_adapt(iexp).av(iav).cycTC,'unif',0)));
    end
end

hm_adapt_norm = cellfun(@(x) x./x(:,1),hm_adapt,'unif',0);
figure
suptitle('SOM cells, first stim responsive')
colormap(brewermap([],'*RdBu'));
for i = 1:2
    subplot(2,2,i)
    [~,sortInd] = sort(hm_adapt_norm{i}(:,8));
    imagesc(hm_adapt_norm{i}(sortInd,:))
    figXAxis([],'Stimulus #',[],1:maxCycles,1:maxCycles)
    figYAxis([],'Cell #',[])
    figAxForm
    title(bxName{i})
    colorbar
%     clim(hm_lim)
end