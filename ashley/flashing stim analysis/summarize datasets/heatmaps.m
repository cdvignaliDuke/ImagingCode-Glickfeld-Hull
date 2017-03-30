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
titleStr = [titleStr mouse(1).expt(1).cells(cellsInd).name];
fnout = fullfile(rc.caOutputDir, dataGroup, [titleStr '_' mouse_str]); 

pre_win = mouse(1).expt(1).win(1).frames;
trans_win = mouse(1).expt(1).win(2).frames;
pre_event_frames = mouse(1).expt(1).pre_event_frames;
post_event_frames = mouse(1).expt(1).post_event_frames;
minTrialLengthFrames = mouse(1).expt(1).info.minTrialLengthFrames;
anti_align = 1;
tar_align = 2;
visual = 1;
auditory = 2;
hits = 1;
cycTime = mouse(1).expt(1).info(1).cyc_time;
cycTimeMs = mouse(1).expt(1).info(1).cyc_time_ms;
nexp = 0;
for imouse = 1:size(mouse,2)
    nexp = nexp+size(mouse(imouse).expt,2);
end
ncells = howManyCells(rc,expt);
ori_stimOn = 11;
%%
doMaxTarResp = 0;
doCellInd = 0;
%% cumulative response for each experment;
trL_cyc = 7;

tc_v_anti = [];
tc_v_anti_err = [];
tc_a_anti = [];
tc_a_anti_err = [];
tc_all_anti = [];

tc_v_tar = [];
tc_v_tar_err = [];

tc_ori = [];
for imouse = 1:size(mouse,2)
    for iexp = 1:size(mouse(imouse).expt,2)
        d_anti = mouse(imouse).expt(iexp).align(anti_align);
        d_tar = mouse(imouse).expt(iexp).align(tar_align);
        
        %anti aligned response
        v = d_anti.av(visual).outcome(hits).cmlvCycResp{trL_cyc};
        a = d_anti.av(auditory).outcome(hits).cmlvCycResp{trL_cyc};
        
        if doCellInd
           ci = mouse(imouse).expt(iexp).cells(cellsInd).ind;
           v = v(:,ci,:);
           a = a(:,ci,:);
        end
        
        tc_v_anti = cat(2,tc_v_anti,mean(v,3));
        tc_v_err = cat(2,tc_v_anti,ste(v,3));
        tc_a_anti = cat(2,tc_a_anti,mean(a,3));
        tc_a_anti_err = cat(2,tc_a_anti,ste(a,3));
        tc_all_anti = cat(2,tc_all_anti,mean(cat(3,v,a),3));
        
        %target aligned response
        
        if doMaxTarResp
            % resp to pref target (abs value)
            v_dir = cellfun(@(x) mean(x,3),d_tar.av(visual).outcome(hits).stimResp,'unif',0);
            [v_dir_max, v_dir_ind] = max(cell2mat(cellfun(@(x) mean(x(trans_win,:),1),v_dir,'unif',0)'),[],1);
            v = nan(size(v_dir{1}));
            nstim = length(v_dir);
            for istim = 1:nstim
                ind = v_dir_ind == istim;
                v_temp = v_dir{istim}(:,ind);
                v(:,ind) = v_temp;            
            end
        else
            v = d_tar.av(visual).outcome(hits).resp;
        end
        
        if doCellInd
            v = v(:,ci,:);
        end
        
        tc_v_tar = cat(2,tc_v_tar,mean(v,3));
        tc_v_tar_err = cat(2,tc_v_tar,ste(v,3));
        
        % tuning aligned response
        
        expt_ind = (strcmp({expt.SubNum},mouse(imouse).expt(iexp).mouse_name)+strcmp({expt.date},mouse(imouse).expt(iexp).date)) == 2;
        dirtuning = expt(expt_ind).dirtuning;
        mName = expt(expt_ind).mouse;
        load(fullfile(rc.ashleyAnalysis,mName,'two-photon imaging',mouse(imouse).expt(iexp).date,dirtuning,'cellsSelect.mat'))
        ori = dFoverF_OriResp_TC;
%         ori = permute(ori,[2,3,1]);
        if doCellInd
            ori = ori(:,ci,:);
        end
        tc_ori = cat(3,tc_ori,ori);
        
    end
end
bl_anti = mean(tc_all_anti(pre_win,:),1);
tc_all_anti_norm = bsxfun(@minus,tc_all_anti,bl_anti);

bl_tar = mean(tc_v_tar(pre_win,:),1);
tc_tar_norm = bsxfun(@minus,tc_v_tar,bl_tar);

tc_ori = tc_ori(:,5:end,:);
%% sort cells according to average response to trial start

oneS_fr = expt(1).frame_rate;

% start_resp = mean(tc_all_anti_norm(trans_win,:),1);
% anti_max_val  = max(tc_all_anti_norm(pre_event_frames:end,:),[],1);
end_resp = mean(tc_all_anti_norm(end-oneS_fr:end,:),1);
% tar_resp = mean(tc_tar_norm(trans_win,:),1);
[resp_sort sort_ind] = sort(end_resp);

% ori_resp = squeeze(mean(tc_ori(:,end-5:end,:),2));
% [max_ori max_ori_ind] = max(ori_resp,[],1);
% sort_ind = [];
% for iori = 1:4
%    ind = find(max_ori_ind == iori);
%    resp_temp = ori_resp(iori,:);
%    
%    [resp_temp_sort sort_temp] = sort(resp_temp(ind));
%    
%    sort_ind = cat(2,sort_ind,ind(sort_temp));
% end



tc_all_sort_anti = tc_all_anti_norm(:,fliplr(sort_ind))';
tc_sort_tar = tc_tar_norm(:,fliplr(sort_ind))';

%% heatmap - anticipation aligned
cb_max = 0.1;
bl_fr = 10;
trL_ind = pre_event_frames - bl_fr:(trL_cyc*cycTime)+pre_event_frames;

tr_tick_fr = (bl_fr:oneS_fr/2:(trL_cyc*cycTime)+bl_fr)+1;
tr_tick_s = 0:0.5:length(tr_tick_fr)-1;

hm=figure;setFigParams4Print('landscape')
suptitle('all trials, all cells, sorted by avg last 1s resp')
colormap(brewermap([],'*RdBu'))

subplot(1,3,1)
anti = imagesc(tc_all_sort_anti(:,trL_ind));
figXAxis(anti.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(anti.Parent);
colorbar
caxis([-cb_max cb_max])
title(['all 100 ms trials >= ' num2str(trL_cyc) ' cycs'])
ylabel('n cells')

%% heatmap - target aligned
max_tar = max(mean(tc_tar_norm(trans_win,:),1));
ms250_fr = ceil(oneS_fr/4);
bl_fr = ms250_fr;
trL_ind_tar = pre_event_frames - bl_fr:bl_fr*2+pre_event_frames;

tr_tick_fr_tar = (0:bl_fr:bl_fr*3)+1;
tr_tick_s_tar = chop(-0.25:0.25:0.5,2);

figure(hm);
subplot(1,3,2)
tar = imagesc(tc_sort_tar(:,trL_ind_tar));
figXAxis(tar.Parent,'time (s)',[],tr_tick_fr_tar,tr_tick_s_tar)
figAxForm(tar.Parent)
colorbar
caxis([-0.2 0.2])
title(['all 100 ms trials, all targets'])
ylabel('n cells')

%% heatmap - orientation tuning

max_ori = max(reshape(squeeze(mean(tc_ori(:,end-5:end,:),2)),size(tc_ori,3)*4,1));
oneS_fr_ori = oneS_fr./10;
stimOnS = round(6./oneS_fr_ori);
trL_ori = size(tc_ori,2);
ori_all = reshape(permute(tc_ori,[2,1,3]),4*trL_ori,[]);
ori_sort = ori_all(:,fliplr(sort_ind))';

tr_block = (0:trL_ori:trL_ori*4);
ori_stim_tick = sort([ ori_stimOn-4:trL_ori:trL_ori*4 tr_block(2:end)]);
ori_stim_label = repmat([0 stimOnS],1,4);

figure(hm);
subplot(1,3,3)
ori = imagesc(ori_sort);
figXAxis(ori.Parent,'time (s)',[],ori_stim_tick,ori_stim_label)
figAxForm(ori.Parent)
hold on
vline(tr_block+0.5,'k-')
colorbar
caxis([-0.8 0.8])
title(['resp drifting gratings; 0, 45, 90, 135, -deg'])
ylabel('n cells')

print([fnout '_heatmap_all_endSort'],'-dpdf','-fillpage');
%% stimuli onset 

%start align
stim1 = 3;
stim_onset = (1:cycTime:cycTime*trL_cyc)+pre_event_frames;
stim_offset = stim_onset+(stim1-1);
stim_on_fr = linspaceNDim(stim_onset,stim_offset,stim1)';

stim_tr = ones(1,pre_event_frames+(cycTime*trL_cyc));
stim_tr(stim_on_fr(:)) = 0;

stim_img = repmat(stim_tr,10,1);


%target align
stim_on_fr = pre_event_frames+(stim1-1);

stim_tr = ones(1,pre_event_frames+(cycTime*trL_cyc));
stim_tr(stim_on_fr(:)) = 0;

stim_img_tar = repmat(stim_tr,10,1);

sm = figure;setFigParams4Print('landscape')
suptitle('stimulus onsets')
colormap gray

%DG
oriStim_fr = 6;
ori_tr = ones(1,size(tc_ori,2));
ori_tr(ori_stimOn-4: ori_stimOn-5+oriStim_fr) = 0;

stim_img_ori = repmat(ori_tr,10,4);

sm = figure;setFigParams4Print('landscape')
suptitle('stimulus onsets')
colormap gray
subplot(1,3,1)
anti = imagesc(stim_img(:,trL_ind));
figXAxis(anti.Parent,'time (s)',[],tr_tick_fr,tr_tick_s);
figAxForm(anti.Parent)
colorbar
title('anti stim')


subplot(1,3,2)
tar = imagesc(stim_img_tar(:,trL_ind_tar));
figXAxis(tar.Parent,'time (s)',[],tr_tick_fr_tar,tr_tick_s_tar);
figAxForm(tar.Parent)
colorbar
title('target stim')

subplot(1,3,3)
ori = imagesc(stim_img_ori);
figXAxis(ori.Parent,'time (s)',[],ori_stim_tick,ori_stim_label);
figAxForm(ori.Parent)
colorbar
title('DG stim')

print([fnout '_heatmap_stim'],'-dpdf','-fillpage');

%% response windows

%response window
resp_tr = ones(1,pre_event_frames+(cycTime*trL_cyc));
resp_tr(trans_win) = 0;

resp_img = repmat(resp_tr,10,1);

rm = figure;setFigParams4Print('landscape')
suptitle('stimulus onsets')
colormap gray

subplot(1,3,1)
anti = imagesc(resp_img(:,trL_ind));
anti.Parent.XTick = tr_tick_fr;
anti.Parent.XTickLabel = tr_tick_s;
colorbar
axis square
title('anti stim')
xlabel('time (s)')
ylabel('n cells')


subplot(1,3,2)
tar = imagesc(resp_img(:,trL_ind_tar));
tar.Parent.XTick = tr_tick_fr_tar;
tar.Parent.XTickLabel = tr_tick_s_tar;
colorbar
axis square
title('target stim')
xlabel('time (s)')
ylabel('n cells')


print([fnout '_heatmap_respWin'],'-dpdf','-fillpage');