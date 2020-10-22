
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
hmLim = [-0.2 0.2];
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
nFr_long = 118;
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
expt_somExptOnly = struct;
expt_sortID = nan(1,nexp);
for im = 1:nMice
    ne = size(mouse(im).expt,2);
    if any(strcmp({som_expt.SubNum},mouse(im).mouse_name))
        for iexp = 1:ne
            if isempty(fieldnames(expt_somExptOnly))
                exptN = 1;
            else
                exptN = exptN+1;
            end
            d = mouse(im).expt(iexp);
            expt_somExptOnly(exptN).expt_name = [mouse(im).mouse_name '-' d.date];
            ind = strcmp({expt.SubNum},mouse(im).mouse_name) & ...
                strcmp({expt.date},d.date);
            if expt(ind).indicator{2}(1:4) == 'flex'
                expt_somExptOnly(exptN).indicator = expt(ind).indicator{2}(6:end);
            else
                expt_somExptOnly(exptN).indicator = expt(ind).indicator{2};
            end
            expt_somExptOnly(exptN).isBx = logical(expt(ind).isBehav);
            expt_sortID(exptN) = find(ind);
            
            if ~isempty(expt(ind).passExpt)
                mouseInd = strcmp({mousePass.mouse_name},mouse(im).mouse_name);
                exptInd = strcmp({mousePass(mouseInd).expt.date},d.date);
                d_pass = mousePass(mouseInd).expt(exptInd);
            else
                d_pass = [];
            end
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
                expt_somExptOnly(exptN).tag(itag).name = tag_name;
                if ~isempty(tag_name)
                    for iav = 1:2
                        if somIsGreen                            
                            expt_somExptOnly(exptN).tag(1).av(iav) = d.tag(2).av(iav);
                            if ~isempty(d_pass)
                                expt_somExptOnly(exptN).tag(1).av_pass(iav) = ...
                                    d_pass.tag(2).av(iav);
                            else
                                expt_somExptOnly(exptN).tag(1).av_pass = [];
                            end
                        else
                            expt_somExptOnly(exptN).tag(itag).av(iav) = d.tag(itag).av(iav);
                            if ~isempty(d_pass)
                                expt_somExptOnly(exptN).tag(itag).av_pass(iav) = ...
                                    d_pass.tag(itag).av(iav);
                                fprintf([num2str(exptN) '\n'])
                            else
                                expt_somExptOnly(exptN).tag(itag).av_pass = [];
                            end
                        end
                    end
                end

            end
        end
    end
end

indicators = unique({expt_somExptOnly.indicator});

% long TC and each stimulus response
cycLengthFr = 11;
resp_SOM = struct;
resp_Other = struct;
for iexp = 1:nexp
    if expt_somExptOnly(iexp).tag(1).name == 'SOM'
        somInd = 1;
        otherInd = 2;
    else
        somInd = 2;
        otherInd = 1;
    end
    % SOM neuron data
    resp_SOM(iexp).isBx = expt_somOnly(iexp).isBx;
    resp_SOM(iexp).channel = channelColor{somInd};
    resp_SOM(iexp).indicator = expt_somOnly(iexp).indicator;
    if ~isempty(expt_somExptOnly(iexp).tag(somInd).av(1).align(1).respTC)
        for iav = 1:2
            d = expt_somExptOnly(iexp).tag(somInd).av(iav).align(1);
            ind_outcome = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
            ind_long = d.nCycles >= max(lateCycles);
            resp_SOM(iexp).av(iav).long = d.respTC(:,:,ind_long) - ...
                mean(d.respTC(basewin,:,ind_long),1);
            cycTC = cell(1,maxCycles);
            lateCycTC = [];
            for icyc = 1:maxCycles
                tc = d.respTC(:,:,d.nCycles >= icyc & ind_outcome);
                cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                cycTC{1,icyc} = tc(...
                    (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                if icyc >= min(lateCycles)
                    lateCycTC = cat(3,lateCycTC,cycTC{1,icyc});
                end
            end
            resp_SOM(iexp).av(iav).cycTC = cycTC;
            resp_SOM(iexp).av(iav).lateCycTC = lateCycTC;
            if isfield(expt_somExptOnly(iexp).tag(somInd),'av_pass')
                fprintf('.\n')
                if ~isempty(expt_somExptOnly(iexp).tag(somInd).av_pass)
                    d = expt_somExptOnly(iexp).tag(somInd).av_pass(iav).align(1);
                    ind_long = d.nCycles >= max(lateCycles);
                    resp_SOM(iexp).av_pass(iav).long = d.respTC(:,:,ind_long) - ...
                        mean(d.respTC(basewin,:,ind_long),1);
                end
            end
        end
    end
    % Other neuron data
    resp_Other(iexp).isBx = expt_somExptOnly(iexp).isBx;
    resp_Other(iexp).channel = channelColor{otherInd};
    resp_Other(iexp).indicator = expt_somExptOnly(iexp).indicator;
    if ~isempty(expt_somExptOnly(iexp).tag(otherInd).av)
        for iav = 1:2
            d = expt_somExptOnly(iexp).tag(otherInd).av(iav).align(1);
            ind_outcome = strcmp(d.outcome,'success')|strcmp(d.outcome,'ignore');
            ind_long = d.nCycles >= max(lateCycles);
            resp_Other(iexp).av(iav).long = d.respTC(:,:,ind_long) - ...
                mean(d.respTC(basewin,:,ind_long),1);
            cycTC = cell(1,maxCycles);
            lateCycTC = [];
            for icyc = 1:maxCycles
                tc = d.respTC(:,:,d.nCycles >= icyc & ind_outcome);
                cycStartOffset = ((icyc-1).*cycLengthFr)+nBaselineFr;
                cycTC{1,icyc} = tc(...
                    (cycStartOffset-nBaselineFr+1):(cycStartOffset+nFrames1s),:,:);
                if icyc >= min(lateCycles)
                    lateCycTC = cat(3,lateCycTC,cycTC{1,icyc});
                end
            end
            resp_Other(iexp).av(iav).cycTC = cycTC;
            resp_Other(iexp).av(iav).lateCycTC = lateCycTC;
            if isfield(expt_somExptOnly(iexp).tag(otherInd),'av_pass')
                if ~isempty(expt_somExptOnly(iexp).tag(otherInd).av_pass)
                    d = expt_somExptOnly(iexp).tag(otherInd).av_pass(iav).align(1);
                    ind_long = d.nCycles >= max(lateCycles);
                    resp_Other(iexp).av_pass(iav).long = d.respTC(:,:,ind_long) - ...
                        mean(d.respTC(basewin,:,ind_long),1);
                end
            end
        end
    end
end

cellInfo_SOM = struct;
cellInfo_Other = struct;
for iexp = 1:nexp
    d_av_long = cat(3,resp_SOM(iexp).av(visualTrials).long,...
        resp_SOM(iexp).av(auditoryTrials).long);
    d_av_first = cat(3,resp_SOM(iexp).av(visualTrials).cycTC{1},...
        resp_SOM(iexp).av(auditoryTrials).cycTC{1});
    d_av_late = cat(3,resp_SOM(iexp).av(visualTrials).lateCycTC,...
        resp_SOM(iexp).av(auditoryTrials).lateCycTC);
    cellInfo_SOM(iexp).responsive_first = ttest(squeeze(mean(d_av_first(respwin,:,:))),...
        squeeze(mean(d_av_first(basewin,:,:))),'dim',2,'tail','right');
    cellInfo_SOM(iexp).responsive_late = ttest(squeeze(mean(d_av_late(respwin,:,:))),...
        squeeze(mean(d_av_late(basewin,:,:))),'dim',2,'tail','right');
    cellInfo_SOM(iexp).responsive_lateWindow = ttest(squeeze(mean(...
        d_av_long((cycLengthFr*length(lateCycles)+1):end,:,:))),...
        squeeze(mean(d_av_long(basewin,:,:))),'dim',2,'tail','right');
    cellInfo_SOM(iexp).suppressed_lateWindow = ttest(squeeze(mean(...
        d_av_long((cycLengthFr*length(lateCycles)+1):end,:,:))),...
        squeeze(mean(d_av_long(basewin,:,:))),'dim',2,'tail','left');
    if ~isempty(resp_Other(iexp).av)
        d_av_long = cat(3,resp_Other(iexp).av(visualTrials).long,...
            resp_Other(iexp).av(auditoryTrials).long);
        d_av_first = cat(3,resp_Other(iexp).av(visualTrials).cycTC{1},...
            resp_Other(iexp).av(auditoryTrials).cycTC{1});
        d_av_late = cat(3,resp_Other(iexp).av(visualTrials).lateCycTC,...
            resp_Other(iexp).av(auditoryTrials).lateCycTC);
        cellInfo_Other(iexp).responsive_first = ttest(squeeze(mean(d_av_first(respwin,:,:))),...
            squeeze(mean(d_av_first(basewin,:,:))),'dim',2,'tail','right');
        cellInfo_Other(iexp).responsive_late = ttest(squeeze(mean(d_av_late(respwin,:,:))),...
            squeeze(mean(d_av_late(basewin,:,:))),'dim',2,'tail','right');
        cellInfo_Other(iexp).responsive_lateWindow = ttest(squeeze(mean(...
            d_av_long((cycLengthFr*length(lateCycles)+1):end,:,:))),...
            squeeze(mean(d_av_long(basewin,:,:))),'dim',2,'tail','right');
        cellInfo_Other(iexp).suppressed_lateWindow = ttest(squeeze(mean(...
            d_av_long((cycLengthFr*length(lateCycles)+1):end,:,:))),...
            squeeze(mean(d_av_long(basewin,:,:))),'dim',2,'tail','left');  
    else
        cellInfo_Other(iexp).responsive_first = [];
    end
end

%% plot heatmaps for naive and training data.
bxInd = cell2mat({resp_SOM.isBx}) == 1;

figure
colormap(brewermap([],'*RdBu'));

d_bx = [];
d_pass = [];
d_naive = [];
for iexp = 1:nexp
    if bxInd(iexp)
        ind = cellInfo_SOM(iexp).responsive_first | ...
            cellInfo_SOM(iexp).responsive_late | ...
            cellInfo_SOM(iexp).responsive_lateWindow | ...
            cellInfo_SOM(iexp).suppressed_lateWindow;
        d_bx = cat(2,d_bx,mean(resp_SOM(iexp).av(visualTrials).long(:,ind,:),3));
        if isfield(resp_SOM(iexp),'av_pass')
            if ~isempty(resp_SOM(iexp).av_pass)
                d_pass = cat(2,d_pass,...
                    mean(resp_SOM(iexp).av_pass(visualTrials).long(:,ind,:),3));
            end
        end
    else
        ind = cellInfo_SOM(iexp).responsive_first | ...
            cellInfo_SOM(iexp).responsive_late | ...
            cellInfo_SOM(iexp).responsive_lateWindow | ...
            cellInfo_SOM(iexp).suppressed_lateWindow;
        d_naive = cat(2,d_naive,mean(resp_SOM(iexp).av(visualTrials).long(:,ind,:),3));
    end
end

subplot 231
[~,sortInd] = sort(mean(d_bx(lateWinFr,:),1));
hm = flipud(d_bx(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('SOM+ - Behav, n=%s',num2str(size(hm,1))))

subplot 232
[~,sortInd] = sort(mean(d_pass(lateWinFr,:),1));
hm = flipud(d_pass(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('SOM+ - Passive, n=%s',num2str(size(hm,1))))

subplot 233
[~,sortInd] = sort(mean(d_naive(lateWinFr,:),1));
hm = flipud(d_naive(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('SOM+ - Naive, n=%s',num2str(size(hm,1))))

d_bx = [];
d_naive = [];
d_pass = [];
for iexp = 1:nexp
    if bxInd(iexp)
        ind = cellInfo_Other(iexp).responsive_first | ...
            cellInfo_Other(iexp).responsive_late | ...
            cellInfo_Other(iexp).responsive_lateWindow | ...
            cellInfo_Other(iexp).suppressed_lateWindow;
        if ~isempty(resp_Other(iexp).av)
            d_bx = cat(2,d_bx,mean(resp_Other(iexp).av(visualTrials).long(:,ind,:),3));
            if isfield(resp_Other(iexp),'av_pass')
                if ~isempty(resp_Other(iexp).av_pass)
                    d_pass = cat(2,d_pass,...
                        mean(resp_Other(iexp).av_pass(visualTrials).long(:,ind,:),3));
                end
            end
        end
    else
        ind = cellInfo_Other(iexp).responsive_first | ...
            cellInfo_Other(iexp).responsive_late | ...
            cellInfo_Other(iexp).responsive_lateWindow | ...
            cellInfo_Other(iexp).suppressed_lateWindow;
        if ~isempty(resp_Other(iexp).av)
            d_naive = cat(2,d_naive,...
                mean(resp_Other(iexp).av(visualTrials).long(:,ind,:),3));
        end
    end
end

subplot 234
[~,sortInd] = sort(mean(d_bx(lateWinFr,:),1));
hm = flipud(d_bx(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('Other - Behav, n=%s',num2str(size(hm,1))))

subplot 235
[~,sortInd] = sort(mean(d_pass(lateWinFr,:),1));
hm = flipud(d_pass(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('Other - Behav, n=%s',num2str(size(hm,1))))

subplot 236
[~,sortInd] = sort(mean(d_naive(lateWinFr,:),1));
hm = flipud(d_naive(:,sortInd)');
imagesc(hm(:,tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr)))
hold on
figXAxis([],'Time from Start (ms)',[],ttLabelFr_long,ttLabel_long)
figYAxis([],'Cell #',[])
figAxForm
colorbar
caxis(hmLim)
title(sprintf('Other - Behav, n=%s',num2str(size(hm,1))))

print([fnout 'suppressedHeatmaps'],'-dpdf','-fillpage')
%% plot TC of suppressed cells
bxInd = cell2mat({resp_SOM.isBx}) == 1;
engage_label = {'Behav','Passive','Naive'};
engage_colors = {[0 0.7 0],[0 0 1],[0 0 0]};
tcSOM = cell(1,3);
tcOther = cell(1,3);
nSOM = zeros(1,3);
nOther = zeros(1,3);
for iexp = 1:nexp
    indSOM = cellInfo_SOM(iexp).suppressed_lateWindow == 1;
    indOther = cellInfo_Other(iexp).suppressed_lateWindow == 1;
    if bxInd(iexp)
        tcSOM{1} = cat(2,tcSOM{1},...
            mean(resp_SOM(iexp).av(visualTrials).long(:,indSOM,:),3));
        nSOM(1) = nSOM(1)+length(indSOM);
        if ~isempty(resp_Other(iexp).av)
            tcOther{1} = cat(2,tcOther{1},...
                mean(resp_Other(iexp).av(visualTrials).long(:,indOther,:),3));
            nOther(1) = nOther(1)+length(indOther);
        end
        if isfield(resp_SOM(iexp),'av_pass')
            if ~isempty(resp_SOM(iexp).av_pass)
                tcSOM{2} = cat(2,tcSOM{2},...
                    mean(resp_SOM(iexp).av_pass(visualTrials).long(:,indSOM,:),3));
                nSOM(2) = nSOM(2)+length(indSOM);
                tcOther{2} = cat(2,tcOther{2},...
                    mean(resp_Other(iexp).av_pass(visualTrials).long(:,indOther,:),3));
                nOther(2) = nOther(2)+length(indOther);
            end
        end
    else
        tcSOM{3} = cat(2,tcSOM{3},...
            mean(resp_SOM(iexp).av(visualTrials).long(:,indSOM,:),3));
        nSOM(3) = nSOM(3)+length(indSOM);        
        if ~isempty(resp_Other(iexp).av)
            tcOther{3} = cat(2,tcOther{3},...
                mean(resp_Other(iexp).av(visualTrials).long(:,indOther,:),3));
            nOther(3) = nOther(3)+length(indOther);
        end
    end
end

figure
subplot 211
leg_lines = [];
leg_label = cell(1,3);
for i = 1:3
    y = mean(tcSOM{i}(tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
    yerr = ste(tcSOM{i}(tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
    hold on
    h = shadedErrorBar_chooseColor(tt_longTC,y,yerr,engage_colors{i});
    leg_lines(i) = h.mainLine;
    leg_label{i} = sprintf('%s: n=%s/%s',engage_label{i},...
        num2str(size(tcSOM{i},2)),num2str(nSOM(i)));
end
legend(leg_lines,leg_label,'location','southwest')
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',[])
figAxForm([],0)
title('SOM+')
hline(0,'k:')

subplot 212
leg_lines = [];
leg_label = cell(1,3);
for i = 1:3
    y = mean(tcOther{i}(tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
    yerr = ste(tcOther{i}(tcStartFrame:((lateCycles(end).*cycLengthFr)+nBaselineFr),:),2);
    hold on
    h = shadedErrorBar_chooseColor(tt_longTC,y,yerr,engage_colors{i});
    leg_lines(i) = h.mainLine;
    leg_label{i} = sprintf('%s: n=%s/%s',engage_label{i},...
        num2str(size(tcOther{i},2)),num2str(nSOM(i)));
end
legend(leg_lines,leg_label,'location','southwest')
figXAxis([],'Time (ms)',[tt_longTC(1) tt_longTC(end)],ttLabel_long,ttLabel_long)
figYAxis([],'dF/F',[])
figAxForm([],0)
title('Other')
hline(0,'k:')

print([fnout 'suppressedTC'],'-dpdf','-fillpage')
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
    early_tc = cat(2,early_tc,mean(resp_SOM(iexp).av(iav).first,3)-...
        mean(mean(resp_SOM(iexp).av(iav).first(basewin_0,:,:),3)));
    late_tc = cat(2,late_tc,mean(resp_SOM(iexp).av(iav).late,3)-...
        mean(mean(resp_SOM(iexp).av(iav).late(basewin_0,:,:),3)));
    ind = cat(1,ind,cellInfo_SOM(iexp).responsive_first);
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
    if resp_SOM(iexp).isBx == 1
        bxInd = 1;
    else
        bxInd = 2;
    end
    early_tc{bxInd} = cat(2,early_tc{bxInd},mean(resp_SOM(iexp).av(iav).first,3)-...
        mean(mean(resp_SOM(iexp).av(iav).first(basewin_0,:,:),3)));
    late_tc{bxInd} = cat(2,late_tc{bxInd},mean(resp_SOM(iexp).av(iav).late,3)-...
        mean(mean(resp_SOM(iexp).av(iav).late(basewin_0,:,:),3)));
    ind{bxInd} = cat(1,ind{bxInd},cellInfo_SOM(iexp).responsive_first);
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
    if strcmp(resp_SOM(iexp).indicator,'GCaMP6s')
        gcampInd = 1;
    else
        gcampInd = 2;
    end
    early_tc{gcampInd} = cat(2,early_tc{gcampInd},mean(resp_SOM(iexp).av(iav).first,3)-...
        mean(mean(resp_SOM(iexp).av(iav).first(basewin_0,:,:),3)));
    late_tc{gcampInd} = cat(2,late_tc{gcampInd},mean(resp_SOM(iexp).av(iav).late,3)-...
        mean(mean(resp_SOM(iexp).av(iav).late(basewin_0,:,:),3)));
    ind{gcampInd} = cat(1,ind{gcampInd},cellInfo_SOM(iexp).responsive_first);
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

%% plot each stimulus adaptation
iav=1;
bxName = {'Behav','Naive'};
hm_lim = [-0.15 0.15];
hm_adapt = cell(1,2);
for iexp = 1:nexp
    if stimResp_adapt(iexp).isBx
        hm_adapt{1} = cat(1,hm_adapt{1},cell2mat(cellfun(@(x)...
            squeeze(mean(mean(x(respwin,cellInfo_SOM(iexp).responsive_first==1,:),1)...
            -mean(x(basewin_0,cellInfo_SOM(iexp).responsive_first==1,:),1),3))',...
            stimResp_adapt(iexp).av(iav).cycTC,'unif',0)));
    else
        hm_adapt{2} = cat(1,hm_adapt{2},cell2mat(cellfun(@(x)...
            squeeze(mean(mean(x(respwin,cellInfo_SOM(iexp).responsive_first==1,:),1)...
            -mean(x(basewin_0,cellInfo_SOM(iexp).responsive_first==1,:),1),3))',...
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