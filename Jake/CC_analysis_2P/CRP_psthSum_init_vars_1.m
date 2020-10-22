
%CRP_psthSum_init_vars_1;

%create variable names
all_trials = [];
all_rew = [];
all_omit = [];
all_unexp = [];
all_earlylick_rew = [];
all_latelick_rew = [];
all_earlylick_omit = [];
all_latelick_omit = [];
all_earlytrial_rew = [];
all_latetrial_rew = [];
all_earlytrial_omit = [];
all_latetrial_omit = [];
all_earlytrial_unexp = [];
all_latetrial_unexp = [];
all_earlytrial_rew_df = [];
all_latetrial_rew_df = [];
all_earlytrial_omit_df = [];
all_latetrial_omit_df = [];
all_earlytrial_unexp_df = [];
all_latetrial_unexp_df = [];
all_earlytrial_rew_lick = [];
all_latetrial_rew_lick = [];
all_earlytrial_omit_lick = [];
all_latetrial_omit_lick = [];
all_earlytrial_unexp_lick = [];
all_latetrial_unexp_lick = [];
all_short_omit = [];
all_long_omit = [];
all_short_unexp = [];
all_long_unexp = [];
all_preomit = [];
all_postomit = [];
all_early_rew_time = zeros(1,nexp);
all_late_rew_time = zeros(1,nexp);
all_early_omit_time = zeros(1,nexp);
all_late_omit_time = zeros(1,nexp);
all_rew_df = [];
all_omit_df = [];
all_unexp_df = [];
all_lick_rew = [];
all_lick_omit = [];
all_lick_unexp = [];
all_postrew_lick_rew = [];
all_postrew_lick_omit = [];
all_postrew_lick_unexp = [];
all_early_postrew_lick_rew = [];
all_late_postrew_lick_rew = [];
all_postrew_lick_rew_lick = [];
all_postrew_lick_omit_lick = [];
all_postrew_lick_unexp_lick = [];
all_early_postrew_lick_rew_lick = [];
all_late_postrew_lick_rew_lick = [];
all_precue_burst = [];
all_precue_single = [];
all_precue_burst_df = [];
all_precue_single_df = [];
all_precue_single_expt = [];
all_precue_burst_expt = [];
all_lastLick_preRew_omit = [];
all_firstLick_postRew_omit = [];
all_lastLick_preRew_rew = [];
all_firstLick_postRew_rew = [];
all_lastLick_preRew_unexp = [];
all_firstLick_postRew_unexp = [];
all_firstPostRewLickEvents = [];
all_expt_bin = [];
all_firstPostRewLickEvents_omit = [];
all_expt_omit_bin = [];
all_firstPostRewLickEvents_unexp = [];
all_expt_unexp_bin = [];
mouse_str = [];
all_area_id = [];
expt_areas = zeros(length(area_list),nexp);
expt_rew_peaks = cell(nexp,length(area_list));
expt_omit_peaks = cell(nexp,length(area_list));
expt_unexp_peaks = cell(nexp,length(area_list));
preresp_rew_range = cell(1,length(area_list));
postresp_rew_range = cell(1,length(area_list));
all_lowlick_rew = [];
all_highlick_rew = [];
all_lowlick_prerew = [];
all_highlick_prerew = [];
all_lowlick_postrew = [];
all_highlick_postrew = [];
all_lowlick_omit = [];
all_highlick_omit = [];
all_lowlick_preomit = [];
all_highlick_preomit = [];
all_lowlick_postomit = [];
all_highlick_postomit = [];
expNums(id).numDendrites = nan(nexp,length(area_list));
expNums(id).avgdFoverF_preRew = nan(nexp,length(area_list));
expNums(id).avgdFoverF_postRew = nan(nexp,length(area_list));
expNums(id).mouse_name = cell(nexp,1);
expNums(id).area_list = area_list;

