ds = '_V1';
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
str = unique({expt.SubNum});
mouse_str = ['i' strjoin(str,'_i')];
load(fullfile(rc.caOutputDir,dataGroup,[mouse_str '_CaSummary' ds '.mat']));
load(fullfile(rc.caOutputDir,dataGroup,[titleStr '_' mouse_str '_modCells.mat']));
% titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name];
fnout = fullfile(rc.caOutputDir, dataGroup,[titleStr '_startAlign']); 

pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
% minTrialLengthFrames = mouse(1).expt(1).info.minTrialLengthFrames;
anti_align = 1;
% tar_align = 2;
visual = 1;
auditory = 2;
hits = 1;
FA = 3;
auroc_comp1 = 3; % 1st stim to target stim
auroc_compL = 1; % 1st stim to target stim
% cycTime = mouse(1).expt(1).info(1).cyc_time;
% cycTimeMs = mouse(1).expt(1).info(1).cyc_time_ms;
nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
% ncells = howManyCells(rc,expt);

%% timing

frRateHz = 30;
fr_2s = 2*frRateHz+pre_event_frames;
fr_bins = (0:frRateHz/2:4*frRateHz)+pre_event_frames; % 500ms bins;
% fr_bins = (0:frRateHz:4*frRateHz)+pre_event_frames; % 1s bins;
% fr_bins = (0:frRateHz*2:4*frRateHz)+pre_event_frames; % 2s bins;

sec_bins = (fr_bins-pre_event_frames)./frRateHz; % time value in s for each bin

%% get trial and cell info aligned

vR1 = [];
aR1 = [];
vTC = [];
aTC = [];
vR_norm_bin = cell(1,length(fr_bins)-1);
aR_norm_bin = cell(1,length(fr_bins)-1);
vR_bin = cell(1,length(fr_bins)-1);
aR_bin = cell(1,length(fr_bins)-1);
vPW_bin = cell(1,length(fr_bins)-1);
aPW_bin = cell(1,length(fr_bins)-1);
vPW_norm_bin = cell(1,length(fr_bins)-1);
aPW_norm_bin = cell(1,length(fr_bins)-1);
vTC_bin = cell(1,length(fr_bins)-1);
aTC_bin = cell(1,length(fr_bins)-1); 
nTr = [];
nTr_bin = cell(1,length(fr_bins)-1); 
v_oneCyc = [];
a_oneCyc = [];
all_oneCyc = [];
cyc100 = [];

ori_pref = [];
untuned = [];
tar_resp = [];
base1_resp = [];
base_sust = [];
auc1_inc = [];
auc1_dec = [];
rst1 = [];
aucL_inc = [];
aucL_dec = [];
rstL = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        cycTimeFr = mouse(imouse).expt(iexp).info.cyc_time;
        dv = mouse(imouse).expt(iexp).align(anti_align).av(visual).outcome(hits);
        da = mouse(imouse).expt(iexp).align(anti_align).av(auditory).outcome(hits);
        v_oneCyc_temp = mean(mouse(imouse).expt(iexp).align(anti_align).av(visual).outcome(FA).oneStim,3);
        a_oneCyc_temp = mean(mouse(imouse).expt(iexp).align(anti_align).av(auditory).outcome(FA).oneStim,3);
        all_oneCyc_temp = mean(cat(3,mouse(imouse).expt(iexp).align(anti_align).av(visual).outcome(FA).oneStim,mouse(imouse).expt(iexp).align(anti_align).av(auditory).outcome(FA).oneStim),3);
        
        all_oneCyc = cat(2,all_oneCyc,all_oneCyc_temp);
        v_oneCyc = cat(2,v_oneCyc,v_oneCyc_temp);
        a_oneCyc = cat(2,a_oneCyc,a_oneCyc_temp);
        if cycTimeFr == 11
            cyc100 = cat(2,cyc100,ones(1,size(v_oneCyc_temp,2)));
        else
            cyc100 = cat(2,cyc100,zeros(1,size(v_oneCyc_temp,2)));
        end
        % resp to first stim and mid-trial length;
        vR = dv.resp;
        aR = da.resp;
        ncells = size(vR,2);
        nTr_temp = (size(vR,3)+size(aR,3))/2;
        vCyc = dv.tcyc;
        aCyc = da.tcyc;
        vR1_temp = squeeze(mean(mean(vR(trans_win,:,:),1),3));
        aR1_temp = squeeze(mean(mean(aR(trans_win,:,:),1),3));
        vR1 = cat(2,vR1,vR1_temp);
        aR1 = cat(2,aR1,aR1_temp);
        vTC = cat(2,vTC,mean(vR(1:fr_2s,:,:),3));
        aTC = cat(2,aTC,mean(aR(1:fr_2s,:,:),3));
        nTr = cat(2,nTr,ones(1,ncells)*nTr_temp);
        
        
        % bin by trial length
        vTrLFr = (vCyc*cycTimeFr)+pre_event_frames;
        aTrLFr = (aCyc*cycTimeFr)+pre_event_frames;
        trLs = cellfun(@(x) size(x,1),dv.cmlvCycResp);
        [~,~,trLs_bin] = histcounts(trLs,fr_bins);      
        for ibin = 1:length(fr_bins)-1;            
            ind = find(trLs_bin == ibin,1,'last');
            if ~isempty(ind)
                vRbin = mean(dv.cmlvCycResp{ind},3);
                aRbin = mean(da.cmlvCycResp{ind},3);
                nTr_temp = (size(dv.cmlvCycResp{ind},3)+size(da.cmlvCycResp{ind},3))/2;
                nTr_bin{ibin} = cat(2,nTr_bin{ibin},ones(1,ncells)*nTr_temp);
                pw = pre_win+((ind-1)*cycTimeFr); %adjust analysis window for tr length
                pwSubCyc = pre_win+((ind-2)*cycTimeFr); % baseline window for previous stim
                tw = trans_win+((ind-1)*cycTimeFr);
                
                vR_norm_temp = mean(vRbin(tw,:),1) - mean(vRbin(pw,:),1);
                aR_norm_temp = mean(aRbin(tw,:),1) - mean(aRbin(pw,:),1);
                vR_temp = mean(vRbin(tw,:),1);
                aR_temp = mean(aRbin(tw,:),1);
                vPW_temp = mean(vRbin(pw,:),1);
                aPW_temp = mean(aRbin(pw,:),1);
                vPW_norm_temp = mean(vRbin(pw,:),1) - mean(vRbin(pwSubCyc,:),1);
                aPW_norm_temp = mean(aRbin(pw,:),1) - mean(aRbin(pwSubCyc,:),1);
                vR_norm_bin{ibin} = cat(2,vR_norm_bin{ibin},vR_norm_temp);
                aR_norm_bin{ibin} = cat(2,aR_norm_bin{ibin},aR_norm_temp);
                vR_bin{ibin} = cat(2,vR_bin{ibin},vR_temp);
                aR_bin{ibin} = cat(2,aR_bin{ibin},aR_temp);
                vPW_norm_bin{ibin} = cat(2,vPW_bin{ibin},vPW_norm_temp);
                aPW_norm_bin{ibin} = cat(2,aPW_bin{ibin},aPW_norm_temp);
                vPW_bin{ibin} = cat(2,vPW_bin{ibin},vPW_temp);
                aPW_bin{ibin} = cat(2,aPW_bin{ibin},aPW_temp);
                        
                ind = vTrLFr >= fr_bins(ibin+1);
                vTC_temp = mean(vR(1:fr_bins(ibin+1),:,ind),3);
                ind = aTrLFr >= fr_bins(ibin+1);
                aTC_temp = mean(aR(1:fr_bins(ibin+1),:,ind),3);

                vTC_bin{ibin} = cat(2,vTC_bin{ibin},vTC_temp);
                aTC_bin{ibin} = cat(2,aTC_bin{ibin},aTC_temp);
            else
                r_nan = nan(1,ncells);
                tc_nan = nan(fr_bins(ibin+1),ncells);
                vR_norm_bin{ibin} = cat(2,vR_norm_bin{ibin},r_nan);
                aR_norm_bin{ibin} = cat(2,aR_norm_bin{ibin},r_nan);
                vR_bin{ibin} = cat(2,vR_bin{ibin},r_nan);
                aR_bin{ibin} = cat(2,aR_bin{ibin},r_nan);
                vTC_bin{ibin} = cat(2,vTC_bin{ibin},tc_nan);
                aTC_bin{ibin} = cat(2,aTC_bin{ibin},tc_nan); 
                nTr_bin{ibin} = cat(2,nTr_bin{ibin},r_nan);
                vPW_bin{ibin} = cat(2,vPW_bin{ibin},r_nan);
                aPW_bin{ibin} = cat(2,aPW_bin{ibin},r_nan);
                vPW_norm_bin{ibin} = cat(2,vPW_bin{ibin},r_nan);
                aPW_norm_bin{ibin} = cat(2,aPW_bin{ibin},r_nan);
            end 
        end
        
        % cell info
%         %tuning
        expt_ind = (strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)+strcmp({expt.date},mouse(imouse).expt(iexp).date)) == 2;
        dirtuning = expt(expt_ind).dirtuning;
        mName = expt(expt_ind).mouse;
        load(fullfile(rc.ashleyAnalysis,mName,'two-photon imaging',mouse(imouse).expt(iexp).date,dirtuning,'cellsSelect.mat'))
        ori_pref = cat(2,ori_pref,ori_ind_all);        
%         %untuned
        temp = zeros(1,ncells);
        temp(mouse(imouse).expt(iexp).cells(6).ind) = 1;
        untuned = cat(2,untuned,temp);
%         %target responsive
        temp = zeros(1,ncells);
        temp(mouse(imouse).expt(iexp).cells(10).ind) = 1;
        tar_resp = cat(2,tar_resp,temp);
%         %first stim responsive
        temp = zeros(1,ncells);
        temp(mouse(imouse).expt(iexp).cells(12).ind) = 1;
        base1_resp = cat(2,base1_resp,temp);
%         %last base responsive
        temp = zeros(1,ncells);
        temp(unique(cat(1,mouse(imouse).expt(iexp).cells(8).ind,mouse(imouse).expt(iexp).cells(9).ind))) = 1;
        base_sust = cat(2,base_sust,temp);
%         %auroc first vs. target
        mi = msModCells(imouse).expt(iexp).av(visual).outcome(hits).mi(1).comp(auroc_comp1).auc(3);
        auc1_inc = cat(2,auc1_inc,mi.value >= 0.5);
        auc1_dec = cat(2,auc1_dec,mi.value < 0.5);
        rst1 = cat(2,rst1,mi.ustat);
%         %auroc last bs vs. target
        mi = msModCells(imouse).expt(iexp).av(visual).outcome(hits).mi(1).comp(auroc_compL).auc(3);
        aucL_inc = cat(2,aucL_inc,mi.value >= 0.5);
        aucL_dec = cat(2,aucL_dec,mi.value < 0.5);
        rstL = cat(2,rstL,mi.rst);

    end
end

%% plotting info
nbins = length(vR_norm_bin);
c_lim = [0 100];
resp_lim_sub = [-0.1 0.25]; 
resp_lim = [-0.3 0.5]; 
resp_lim_norm = [-10 10]; 
np1 = 3;
np2 = 3;
%% plot scatters of respones for each time bin 
cell_ind = logical(base1_resp | base_sust | tar_resp);
% cell_ind = 1:length(vR1);
% PREV RESPONSE SUBTRACTED
figure; setFigParams4Print('landscape');
suptitle('response to last vis stim, previous stim subtracted')
subplot(np1,np2,1)
s = scatter(vR1(cell_ind),aR1(cell_ind),100,nTr(cell_ind),'.');
hold on
plot(resp_lim_sub,resp_lim_sub,'k--')
figXAxis(s.Parent,'visual',resp_lim_sub)
figYAxis(s.Parent,'auditory',resp_lim_sub)
figAxForm(s.Parent)
title('first stim')
colorbar
caxis(c_lim)
for ibin = 1:nbins
    subplot(np1,np2,ibin+1)
    s = scatter(vR_norm_bin{ibin}(cell_ind),aR_norm_bin{ibin}(cell_ind),100,nTr_bin{ibin}(cell_ind),'.');
    hold on
    plot(resp_lim_sub,resp_lim_sub,'k--')
    figXAxis(s.Parent,'visual',resp_lim_sub)
    figYAxis(s.Parent,'auditory',resp_lim_sub)
    figAxForm(s.Parent)
    title(['< ' num2str(sec_bins(ibin+1)) 's'])
    colorbar
    caxis(c_lim)
end
print([fnout '_scat_lastBS2previous'],'-dpdf','-fillpage')

% RELATIVE TO START
figure; setFigParams4Print('landscape');
suptitle('response to last vis stim relative to start')
subplot(np1,np2,1)
s = scatter(vR1(cell_ind),aR1(cell_ind),100,nTr(cell_ind),'.');
hold on
plot(resp_lim,resp_lim,'k--')
figXAxis(s.Parent,'visual',resp_lim)
figYAxis(s.Parent,'auditory',resp_lim)
figAxForm(s.Parent)
title('first stim')
colorbar
caxis(c_lim)
for ibin = 1:nbins
    subplot(np1,np2,ibin+1)
    s = scatter(vR_bin{ibin}(cell_ind),aR_bin{ibin}(cell_ind),100,nTr_bin{ibin}(cell_ind),'.');
    hold on
    plot(resp_lim,resp_lim,'k--')
    figXAxis(s.Parent,'visual',resp_lim)
    figYAxis(s.Parent,'auditory',resp_lim)
    figAxForm(s.Parent)
    title(['< ' num2str(sec_bins(ibin+1)) 's'])
    colorbar
    caxis(c_lim)
end
print([fnout '_scat_lastBS2start'],'-dpdf','-fillpage')

% BASELINE RELATIVE TO START
figure; setFigParams4Print('landscape');
suptitle('baseline dF/F relative to start')

for ibin = 1:nbins
    subplot(np1,np2,ibin+1)
    s = scatter(vPW_bin{ibin}(cell_ind),aPW_bin{ibin}(cell_ind),100,nTr_bin{ibin}(cell_ind),'.');
    hold on
    plot(resp_lim,resp_lim,'k--')
    figXAxis(s.Parent,'visual',resp_lim)
    figYAxis(s.Parent,'auditory',resp_lim)
    figAxForm(s.Parent)
    title(['< ' num2str(sec_bins(ibin+1)) 's'])
    colorbar
    caxis(c_lim)
end
print([fnout '_scat_preStimWindow'],'-dpdf','-fillpage')

% BASELINE RELATIVE TO PREV STIM
figure; setFigParams4Print('landscape');
suptitle('baseline dF/F relative to prev stim')

for ibin = 1:nbins
    subplot(np1,np2,ibin+1)
    s = scatter(vPW_norm_bin{ibin}(cell_ind),aPW_norm_bin{ibin}(cell_ind),100,nTr_bin{ibin}(cell_ind),'.');
    hold on
    plot(resp_lim,resp_lim,'k--')
    figXAxis(s.Parent,'visual',resp_lim)
    figYAxis(s.Parent,'auditory',resp_lim)
    figAxForm(s.Parent)
    title(['< ' num2str(sec_bins(ibin+1)) 's'])
    colorbar
    caxis(c_lim)
end
print([fnout '_scat_preStimWindowNorm'],'-dpdf','-fillpage')
%% heatmaps of different trial lengths
% fig params
cb_max = 0.2;
bl_fr = 10;
tr_tick_fr = bl_fr:frRateHz/2:4*frRateHz;
tr_tick_s = (0:frRateHz/2:4*frRateHz)/frRateHz;
tr_start = pre_event_frames-bl_fr;

% sorting
% allTC = vTC-aTC;
% meanSub = mean(allTC(pre_event_frames:end,:),1);
% [~,sort_ind] = sort(meanSub);
allTC = mean(cat(3,vTC,aTC),3);
maxTC = max(allTC(pre_event_frames:end,:),[],1);
[~,sort_ind] = sort(maxTC);
sort_ind = fliplr(sort_ind);


v_hm = figure; setFigParams4Print('landscape'); suptitle({'vis trials';'sorted by avg resp to first stim'})
colormap(brewermap([],'*RdBu'))
a_hm = figure; setFigParams4Print('landscape'); suptitle({'aud trials';'sorted by avg resp to first stim'})
colormap(brewermap([],'*RdBu'))
sub_hm = figure; setFigParams4Print('landscape'); suptitle({'vis - aud';'sorted by avg resp to first stim'})
colormap(brewermap([],'*RdBu'))

for ibin = 1:nbins
    vTC_sort = vTC_bin{ibin}(tr_start:end,sort_ind)';
    trL = size(vTC_sort,2);
    tick_ind = tr_tick_fr <= trL;
    tick_fr  = tr_tick_fr(tick_ind);
    tick_s = tr_tick_s(tick_ind);
    cell_ind = ~isnan(mean(vTC_sort,2));

    figure(v_hm);
    subplot(np1,np2,ibin+1)
    hm = imagesc(vTC_sort(cell_ind,:));  
    figXAxis(hm.Parent,'time (s)',[],tick_fr,tick_s);
    figAxForm(hm.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['< ' num2str(sec_bins(ibin+1)) 's'])
    ylabel('n cells')

    aTC_sort = aTC_bin{ibin}(tr_start:end,sort_ind)';

    figure(a_hm);
    subplot(np1,np2,ibin+1)
    hm = imagesc(aTC_sort(cell_ind,:));  
    figXAxis(hm.Parent,'time (s)',[],tick_fr,tick_s);
    figAxForm(hm.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['< ' num2str(sec_bins(ibin+1)) 's'])
    ylabel('n cells')

    subTC_sort = vTC_sort(cell_ind,:)-aTC_sort(cell_ind,:);

    figure(sub_hm);
    subplot(np1,np2,ibin+1)
    hm = imagesc(subTC_sort);  
    figXAxis(hm.Parent,'time (s)',[],tick_fr,tick_s);
    figAxForm(hm.Parent);
    colorbar
    caxis([-cb_max cb_max])
    title(['< ' num2str(sec_bins(ibin+1)) 's'])
    ylabel('n cells')
end        
        