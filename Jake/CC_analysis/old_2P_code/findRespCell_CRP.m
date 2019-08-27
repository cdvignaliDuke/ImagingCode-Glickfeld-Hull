function [NR_h, NR_p, NR_resp_cells, NR_resp_avg, NR_resp_sem, NR_base, NR_resp] = findRespCell_CRP(NR_movie, ops)
%function for identifying significantly responsive neurons using the df/f traces from dendritic masks
%NR_movie   dim1=trials  dim2=cells  dim3=frames
% 
% ops.base_cue_buffer 
% ops.resp_cue_buffer 
% ops.effect_sign 
% ops.event_type
% ops.ifi
% ops.cue_frame
% ops.rew_frame

if isempty(NR_movie)
    NR_h =[]; NR_p =[]; NR_resp_cells =[]; NR_resp_avg =[]; NR_resp_sem =[]; NR_base =[]; NR_resp =[];
    return
end
ifi = ops.ifi;
effect_sign = ops.effect_sign;

%define baseline and response windows
base_window = ops.cue_frame-round(ops.base_cue_buffer./double(ifi)):ops.cue_frame-1;  %cue-500ms:cue
%base_window2 = ops.rew_frame-round(ops.base_cue_buffer./double(ifi)) : ops.rew_frame-1; %rew-500ms:rew
base_window2 = ops.rew_frame-round(300./double(ifi)) : ops.rew_frame-1; %rew-300ms:rew
if strcmp(ops.event_type, 'reward')
    resp_window = ops.rew_frame : ops.rew_frame + round(ops.resp_cue_buffer./double(ifi));
elseif strcmp(ops.event_type, 'cue')
    resp_window = ops.cue_frame+double(round(200/ifi)):ops.cue_frame+double(round(200/ifi))+round(ops.resp_cue_buffer./double(ifi));
end

%find the max value indeces in the base/resp windows
[~, ind_resp] = max(nanmean(NR_movie(:,:,resp_window),1),[],3);
[~, ind_resp_neg] = min(nanmean(NR_movie(:,:,resp_window),1),[],3);
[~, ind_base] = max(nanmean(NR_movie(:,:,base_window),1),[],3);
[~, ind_base2] = max(nanmean(NR_movie(:,:,base_window2),1),[],3);

%allocate memory and define variables
nCells = size(NR_movie, 2);
NR_base = zeros(size(NR_movie,1),nCells);
NR_base2 = zeros(size(NR_movie,1),nCells);
NR_resp = zeros(size(NR_movie,1),nCells);
NR_resp_neg = zeros(size(NR_movie,1),nCells);

%for each cell ID both baselines, pos resp, and negative resp
for ic = 1:nCells
    %find mean baseline response
%     NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,[base_window]),3)); 
%     NR_base2(:,ic) = squeeze(mean(NR_movie(:,ic,[base_window2]),3));  
    %define baseline window and store value 
    if base_window(ind_base(ic))+round(50./double(ifi)) > size(base_window,2) %make sure there is enough room in the rest of the TC for the window. 
        NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_window(ind_base(ic))-round(100./double(ifi)):base_window(ind_base(ic))),3));
    else
        NR_base(:,ic) = squeeze(mean(NR_movie(:,ic,base_window(ind_base(ic))-round(50./double(ifi)):base_window(ind_base(ic))+round(50./double(ifi))),3));
    end
    %define baseline2 window and store value 
    NR_base2(:,ic) = squeeze(mean(NR_movie(:,ic,base_window2),3));
%     if base_window2(ind_base2(ic))+round(50./double(ifi)) > size(base_window2,2) %make sure there is enough room in the rest of the TC for the window. 
%         NR_base2(:,ic) = squeeze(mean(NR_movie(:,ic,base_window2(ind_base2(ic))-round(100./double(ifi)):base_window2(ind_base2(ic))),3));
%     else
%         NR_base2(:,ic) = squeeze(mean(NR_movie(:,ic,base_window2(ind_base2(ic))-round(50./double(ifi)):base_window2(ind_base2(ic))+round(50./double(ifi))),3));
%     end
    
    %define response window and store value  POSITIVE
    %     NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,[resp_window]),3));
    if strcmp(ops.event_type, 'reward')
        NR_resp(:,ic) = squeeze(NR_movie(:,ic,resp_window(ind_resp(ic))));
    elseif resp_window(ind_resp(ic))+round(50./double(ifi)) > size(resp_window,2) %make sure there is enough room in the rest of the TC for the window. 
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_window(ind_resp(ic))-round(100./double(ifi)):resp_window(ind_resp(ic))),3));
    else
        NR_resp(:,ic) = squeeze(mean(NR_movie(:,ic,resp_window(ind_resp(ic))-round(50./double(ifi)):resp_window(ind_resp(ic))+round(50./double(ifi))),3));
    end
    
    %define response window and store value  NEGATIVE
    %     NR_resp_neg(:,ic) = squeeze(mean(NR_movie(:,ic,[resp_window]),3));
    if resp_window(ind_resp_neg(ic))+round(50./double(ifi)) > size(resp_window,2) %make sure there is enough room in the rest of the TC for the window. 
        NR_resp_neg(:,ic) = squeeze(mean(NR_movie(:,ic,resp_window(ind_resp_neg(ic))-round(100./double(ifi)):resp_window(ind_resp_neg(ic))),3));
    else
        NR_resp_neg(:,ic) = squeeze(mean(NR_movie(:,ic,resp_window(ind_resp_neg(ic))-round(50./double(ifi)):resp_window(ind_resp_neg(ic))+round(50./double(ifi))),3));
    end
end

%% 2. ttest for significant responses
%if all values are nans then h and p are zeroes - could happen if you dont have enough trials? 
if isempty(find(~isnan(NR_base)))
    NR_h = zeros(1,size(NR_base,2));
    NR_p = zeros(1,size(NR_base,2));
else
    if strcmp(effect_sign, 'pos') ==1
        ttest_dir = 'right';
        if strcmp(ops.event_type, 'reward')
            %[NR_h, NR_p] = ttest(NR_resp, NR_base2, 'dim', 1, 'tail', 'right', 'alpha', 0.05);
            [NR_h, NR_p] = ttest_forloop(NR_resp, NR_base2, ttest_dir);
        elseif strcmp(ops.event_type, 'cue')
            %[NR_h, NR_p] = ttest(NR_resp, NR_base, 'dim', 1, 'tail', 'right', 'alpha', 0.05);
            [NR_h, NR_p] = ttest_forloop(NR_resp, NR_base, ttest_dir);
        end
    elseif strcmp(effect_sign, 'neg') ==1
        ttest_dir = 'left';
        if strcmp(ops.event_type, 'reward')
            %[NR_h1, NR_p1] = ttest(NR_resp_neg, NR_base2, 'dim', 1, 'tail', 'left', 'alpha', 0.05);
            [NR_h1, NR_p1] = ttest_forloop(NR_resp_neg, NR_base2, ttest_dir);
            %[NR_h2, NR_p2] = ttest(NR_resp_neg, NR_base, 'dim', 1, 'tail', 'left', 'alpha', 0.05);
            [NR_h2, NR_p2] = ttest_forloop(NR_resp_neg, NR_base, ttest_dir);
            NR_h1(NR_h1~=1) = 0;  NR_h2(NR_h2~=1) = 0;
            NR_h = NR_h1 & NR_h2;
            NR_p = [NR_p1; NR_p2];
        elseif strcmp(ops.event_type, 'cue')
            [NR_h, NR_p] = ttest(NR_resp_neg, NR_base, 'dim', 1, 'tail', 'left', 'alpha', 0.05);
            [NR_h, NR_p] = ttest_forloop(NR_resp_neg, NR_base, ttest_dir);
        end
    end
end

%% 3 set a minimum magnitude for the difference between base and resp
NR_resp_diff = (nanmean(NR_resp)-nanmean(NR_base)) > 0.015;
NR_resp_diff2 = (nanmean(NR_resp)-nanmean(NR_base2)) > 0.015;
NR_resp_diff(NR_resp_diff~=1) = 0;
NR_resp_diff2(NR_resp_diff2~=1) = 0;
NR_h(NR_h~=1) = 0;
if strcmp(ops.event_type, 'reward')
    NR_h = NR_h & NR_resp_diff2;
elseif strcmp(ops.event_type, 'cue')
    NR_h = NR_h & NR_resp_diff;
end

%% 4. define cells by response to event
NR_h = logical(NR_h);
NR_resp_cells = find(NR_h);

NR_resp_avg = []; %nanmean((NR_resp-NR_base),1);
NR_resp_sem = [];  %nanstd((NR_resp-NR_base),[],1)./sqrt(size(NR_resp,1));

end

function [this_h, this_p] = ttest_forloop(resp_samp, base_samp, ttest_dir)
this_h = [];
this_p = [];
for cell_num = 1:size(resp_samp,2)
    this_cell_resp = resp_samp(:,cell_num);
    this_cell_base = base_samp(:,cell_num);
    [this_h(cell_num), this_p(cell_num)] = ttest(this_cell_resp(isfinite(this_cell_resp)), this_cell_base(isfinite(this_cell_base)), 'dim', 1, 'tail', ttest_dir, 'alpha', 0.05);
end
end



