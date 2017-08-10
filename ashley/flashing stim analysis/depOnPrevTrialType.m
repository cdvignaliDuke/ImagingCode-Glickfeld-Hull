% part of FSAVsum_anti

%% create structure of visual and auditory trial types 
% example: visual trials will be 1x3 cell array that's a cumulative
% response for all trials, prev same trials, and prev different trials,
% corresponding  naming array also made here
minTrNexpt = 2;

v_prev_success = strcmp(dv.prevTrOut,'success') | strcmp(dv.prevTrOut,'ignore');
v_prev2_success = strcmp(dv.prev2TrOut,'success') | strcmp(dv.prev2TrOut,'ignore');
a_prev_success = strcmp(da.prevTrOut,'success') | strcmp(da.prevTrOut,'ignore');
a_prev2_success = strcmp(da.prev2TrOut,'success') | strcmp(da.prev2TrOut,'ignore');

v_prev_ind = cell(1,6);
a_prev_ind = cell(1,6);

% previous trial visual
p_ind = dv.prevTrType == 0 & trl_ind_v & v_prev_success;
v_prev_ind{1} = p_ind;
p_ind = da.prevTrType == 0 & trl_ind_a & a_prev_success;
a_prev_ind{1} = p_ind;

% previous trial auditory
p_ind = dv.prevTrType == 1 & trl_ind_v & v_prev_success;
v_prev_ind{2} = p_ind;
p_ind = da.prevTrType == 1 & trl_ind_a & a_prev_success;
a_prev_ind{2} = p_ind;

% previous 2 trials visual
p_ind = dv.prev2TrType == 0 & dv.prevTrType == 0 & trl_ind_v & v_prev_success & v_prev2_success;
v_prev_ind{3} = p_ind;
p_ind = da.prev2TrType == 0 & da.prevTrType == 0 & trl_ind_a & a_prev_success & a_prev2_success;
a_prev_ind{3} = p_ind;

% previous 2 trial auditory
p_ind = dv.prev2TrType == 1 & dv.prevTrType == 1 & trl_ind_v & v_prev_success & v_prev2_success;
v_prev_ind{4} = p_ind;
p_ind = da.prev2TrType == 1 & da.prevTrType == 1 & trl_ind_a & a_prev_success & a_prev2_success;
a_prev_ind{4} = p_ind;

% previous 2 trial vis->aud
p_ind = dv.prev2TrType == 0 & dv.prevTrType == 1 & trl_ind_v & v_prev_success & v_prev2_success;
v_prev_ind{5} = p_ind;
p_ind = da.prev2TrType == 0 & da.prevTrType == 1 & trl_ind_a & a_prev_success & a_prev2_success;
a_prev_ind{5} = p_ind;

% previous 2 trial aud->vis
p_ind = dv.prev2TrType == 1 & dv.prevTrType == 0 & trl_ind_v & v_prev_success & v_prev2_success;
v_prev_ind{6} = p_ind;
p_ind = da.prev2TrType == 1 & da.prevTrType == 0 & trl_ind_a & a_prev_success & a_prev2_success;
a_prev_ind{6} = p_ind;

trl = fr_bins(ibin+1);
nc = size(dv.resp,2);
for imat = 1:6
    p_ind = v_prev_ind{imat};
    if sum(p_ind) >= minTrNexpt
        vp_temp = mean(dv.resp(1:fr_bins(ibin+1),:,p_ind),3);
    else
        vp_temp = NaN(trl,nc);
    end
    vPrev_mat{imat} = cat(2,vPrev_mat{imat},vp_temp);
    
    p_ind = a_prev_ind{imat};
    if sum(p_ind) >= minTrNexpt
        ap_temp = mean(da.resp(1:fr_bins(ibin+1),:,p_ind),3);
    else        
        ap_temp = NaN(trl,nc);
    end
    aPrev_mat{imat} = cat(2,aPrev_mat{imat},ap_temp);
end

%% repeated trial - mod cells
if strcmp(ds,'_V1')
    v_2rep_early = squeeze(mean(dv.resp(early_win,:,v_prev_ind{1}),1));
    a_2rep_early = squeeze(mean(da.resp(early_win,:,a_prev_ind{2}),1));
    v_2rep_late = squeeze(mean(dv.resp(late_win,:,v_prev_ind{1}),1));
    a_2rep_late = squeeze(mean(da.resp(late_win,:,a_prev_ind{2}),1));

    H_av_early_2rep_temp = ttest2(v_2rep_early, a_2rep_early, 'tail', 'both','dim',2);
    H_av_late_2rep_temp = ttest2(v_2rep_late, a_2rep_late, 'tail', 'both', 'dim', 2);

    H_av_early_2rep = cat(1,H_av_early_2rep,H_av_early_2rep_temp);
    H_av_late_2rep = cat(1,H_av_late_2rep,H_av_late_2rep_temp);
end