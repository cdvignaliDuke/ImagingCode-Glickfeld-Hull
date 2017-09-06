clear all
close all
ds = '';
%%
rc = behavConstsAV;
if strcmp(rc.name,'ashle')
    dataGroup = ['awFSAVdatasets' ds];
else
    dataGroup = [];
end
eval(dataGroup)
titleStr = ds;
if strcmp(titleStr, '')
    titleStr = 'V1_100ms';
else
    titleStr = titleStr(2:end);
end

load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' ds '.mat']));
load(fullfile(rc.caOutputDir,dataGroup,[titleStr '_' mouse_str '_modCells.mat']));
fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_startAlign']); 

pressAlign = 1;
visualTrials = 1;
auditoryTrials = 2;
hitTrials = 1;
missTrials = 2;
cycLengthFr = mouse(1).expt(1).info.cyc_time;
frRateHz = expt(1).frame_rate;
onTimeFr = 0.1*frRateHz;
nMonitorDelayFr = 2;
nBaselineFr = mouse(1).expt(1).pre_event_frames;
timestamp_1cyc = (-nBaselineFr:cycLengthFr-1);
plotTimeRange_1cyc = [-5 chop(timestamp_1cyc(end),2)];
nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end

% response to each of first 5 cyc
basewin = [32 34];
respwin = [36 38];
nCycles = 8;
nav = 2;
tcEaCellEaCyc = cell(nCycles,nexp,nav);
meanRespEaCell = cell(nexp,nav);
nTrialsEaCyc = zeros(nCycles,nexp,nav);
exptName = cell(1,nexp);
h = cell(1,nexp);
p = cell(nexp,nCycles);
good_ind = cell(nexp,nCycles);
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
            if imouse == 1 && iexp == 1
                exptN = 1;
            else
                exptN = exptN+1;
            end
        exptName{exptN} = [mouse(imouse).expt(iexp).mouse_name '-' ...
            mouse(imouse).expt(iexp).date];
        for iav = 1:nav
            d = mouse(imouse).expt(iexp).align(pressAlign).av(iav);        
            for icyc = 1:nCycles
                tcThisCyc = cat(3,d.outcome(hitTrials).cmlvCycResp{icyc}...
                    (:,:,:),d.outcome(missTrials).cmlvCycResp...
                    {icyc}(:,:,:));
                if iav == 1
                    if icyc == 1
                        base_resp = mean(tcThisCyc(basewin(1):basewin(2),:,:),1);
                    end
                    good_ind{exptN,icyc} = zeros(1,size(tcThisCyc,2));
                    h{exptN} = zeros(1,size(tcThisCyc,2));
                    p{exptN,icyc} = zeros(1,size(tcThisCyc,2));
                    for iCell = 1:size(tcThisCyc,2)
                        if icyc == 1
                            [h{exptN}(1,iCell) p{exptN,icyc}(1,iCell)] = ttest(squeeze(mean(tcThisCyc(respwin(1):respwin(2),iCell,:),1)), squeeze(base_resp(:,iCell,:)),'tail','right');
                            if h{exptN}(1,iCell) == 1
                                [max_val max_ind{exptN}(1,iCell)] = max(diff(mean(tcThisCyc(26:end,iCell,:),3),1),[],1);
                                if max_ind{exptN}(1,iCell) < respwin(end)-25
                                    good_ind{exptN,icyc}(1,iCell) = 1;
                                end
                            end
                        else
                            [good_ind{exptN,icyc}(1,iCell) p{exptN,icyc}(1,iCell)] = ttest2(squeeze(mean(tcThisCyc(respwin(1):respwin(2),iCell,:),1)), squeeze(base_resp(:,iCell,:)),'tail','right');
                        end
                    end
                end
                thisExptRespEaCyc = nan(nCycles,size(tcThisCyc,2));
                extraFrames = cycLengthFr*(icyc-1);
                responseAllTrials = mean(tcThisCyc(...
                    responseWindow+extraFrames,:,:),1);
                baselineAllTrials = mean(tcThisCyc(...
                    baselineWindow+extraFrames,:,:),1);
                thisExptRespEaCyc(icyc,:) = squeeze(mean(responseAllTrials - ...
                    baselineAllTrials,3));
                nTrialsEaCyc(icyc,exptN,iav) = size(responseAllTrials,3);
                tcEaCellEaCyc{icyc,exptN,iav} = mean(bsxfun(@minus,tcThisCyc(...
                    end-(cycLengthFr+nBaselineFr)+1:end,:,:),...
                    baselineAllTrials),3);

            end
            meanRespEaCell{exptN,iav} = thisExptRespEaCyc;
        end
    end
end
%% figures

figure; 
for i = 1:nexp
    subplot(4,5,i);
    ind = find(good_ind{i,1});
    nCells = length(ind);
    for ii = 1:5
        plot(mean(tcEaCellEaCyc{ii,i,1}(25:end,ind),2)); 
        hold on; 
        title(num2str(length(ind)))
    end 
end


figure
subplot(2,2,1)
temp_tcs = cell(nCycles,nav);
for ii = 1:nCycles
    for i = 1:nexp
        for iav = 1:nav
            ind = find(good_ind{i,1});
            temp_tcs{ii,iav} = [temp_tcs{ii,iav} tcEaCellEaCyc{ii,i,iav}(:,ind)]; 
        end
    end
    if ii == 1
        big_ind = find(mean(temp_tcs{1,1}(respwin,:),1)>0.002);
    end
    plot(((26:size(temp_tcs{ii,1},1))-33).*(1000/frRateHz),mean(temp_tcs{ii,1}(26:end,big_ind),2))
    hold on
    for iav = 1:nav
        avg_resp(ii,iav) = mean(mean(temp_tcs{ii,iav}(respwin,big_ind),1),2);
        sem_resp(ii,iav) = ste(mean(temp_tcs{ii,iav}(respwin,big_ind),1),2);
    end
end
title(['n = ' num2str(length(big_ind))])
xlabel('Time (ms)')
ylabel('dF/F')
vline((respwin-33).*(1000/frRateHz))

subplot(2,2,3)
for iav = 1:nav
    errorbar(1:nCycles, avg_resp(:,iav), sem_resp(:,iav));
    hold on
end
xlim([0 nCycles+1])
xlabel('Stimulus number')
ylabel('dF/F')


subplot(2,2,2)
avg_tcs = zeros(size(temp_tcs{1,1},1),nCycles,nav);
norm_tcs = cell(nCycles,nav);
avg_resp_norm = zeros(nCycles,nav);
sem_resp_norm = zeros(nCycles,nav);
for ii = 1:nCycles
    for iav = 1:nav
        norm_tcs{ii,iav} = temp_tcs{ii,iav}(:,big_ind)./mean(temp_tcs{1,1}(respwin,big_ind),1);
        avg_tcs(:,ii,iav) = mean(norm_tcs{ii,iav},2);
        avg_resp_norm(ii,iav) = mean(mean(norm_tcs{ii,iav}(respwin,:),1),2);
        sem_resp_norm(ii,iav) = ste(mean(norm_tcs{ii,iav}(respwin,:),1),2);
        if iav == 1
            plot(((26:size(avg_tcs,1))-33).*(1000/frRateHz),avg_tcs(26:end,ii,iav))
            hold on
        end
    end
end
title(['n = ' num2str(length(big_ind))])
xlabel('Time (ms)')
ylabel('Normalized dF/F')
vline((respwin-33).*(1000/frRateHz))

subplot(2,2,4)
for iav = 1:nav
    errorbar(1:nCycles, avg_resp_norm(:,iav), sem_resp_norm(:,iav));
    hold on
end
ylim([0 1.1])
xlim([0 nCycles+1])
xlabel('Stimulus number')
ylabel('Normalized dF/F')

print('Z:\home\lindsey\Analysis\2P\Adaptation\Adaptation_allCells_AV_AWdataset.pdf','-dpdf','-bestfit')


figure; 
colmat = ['k','c'];
norm_tcs_cum = cell(1,nav);
temp_tcs_cum = cell(1,nav);
for ii = 1:nCycles
    subplot(3,4,ii)
    for iav = 1:nav
        shadedErrorBar(((26:size(avg_tcs,1))-33).*(1000/frRateHz), mean(norm_tcs{ii,iav}(26:end,:),2), ste(norm_tcs{ii,iav}(26:end,:),2),colmat(:,iav));
        hold on
        if ii > 2 & ii < 7
            norm_tcs_cum{iav} = cat(3, norm_tcs_cum{iav}, norm_tcs{ii,iav});
            temp_tcs_cum{iav} = cat(3, temp_tcs_cum{iav}, temp_tcs{ii,iav}(:,big_ind));
        end
    end
    title(['Cycle #' num2str(ii)])
end
subplot(3,2,5)
for iav = 1:nav
    shadedErrorBar(((26:size(avg_tcs,1))-33).*(1000/frRateHz), mean(mean(norm_tcs_cum{iav}(26:end,:,:),3),2), ste(mean(norm_tcs{ii,iav}(26:end,:,:),3),2),colmat(:,iav));
    hold on
end
title(['Cycle #3-6'])

subplot(3,2,6)
scatter(mean(mean(temp_tcs_cum{1}(respwin,:,:),1),3),mean(mean(temp_tcs_cum{2}(respwin,:,:),1),3),'ok')
[h_ppr p_ppr] = ttest(mean(mean(temp_tcs_cum{1}(respwin,:,:),1),3),mean(mean(temp_tcs_cum{2}(respwin,:,:),1),3));
xlabel('Visual')
ylabel('Auditory')
xlim([-0.01 0.05])
ylim([-0.01 0.05])
refline(1,0) 
axis square
title(['p = ' num2str(chop(p_ppr,2))])

print('Z:\home\lindsey\Analysis\2P\Adaptation\Adaptation_allCellsByCycle_AV_AWdataset.pdf','-dpdf','-bestfit')

