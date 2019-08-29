clear all
for id = [1 2 3]
    close all
    CRP_expt_list_Crus2019
    lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
    behav_dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\Behavior\Data';
    jake_dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';
    nexp = size(expt(id).date,1);
    for iexp = 1:nexp
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];
%         if ~exist(fullfile(lg_out,img_fn, [img_fn '_ROI_TCs.mat']))
%             load(fullfile(jake_dir,img_fn, [img_fn '_ROI_TCs.mat']))
%             save(fullfile(lg_out,img_fn, [img_fn '_ROI_TCs.mat']), 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
%         end
%         if ~exist(fullfile(lg_out,img_fn, [img_fn '_input.mat']))
%             load(fullfile(jake_dir,img_fn, [img_fn '_input.mat']))
%             save(fullfile(lg_out,img_fn, [img_fn '_input.mat']), 'input');
%         end
%         if ~exist(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))
%             load(fullfile(jake_dir,img_fn, [img_fn '_targetAlign.mat']))
%             save(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']), 'ind_rew', 'ind_omit', 'ind_unexp', 'targetAlign_events', 'targetAligndFoverF', 'prewin_frames', 'postwin_frames', 'tt', 'frameRateHz', 'ind_omit_short','ind_omit_long','ind_unexp_short','ind_unexp_long','ind_rew_preomit','ind_rew_postomit');
%         end
%         if ~exist(fullfile(lg_out,img_fn, [img_fn '_cueAlignLick.mat']))
%             load(fullfile(jake_dir,img_fn, [img_fn '_cueAlignLick.mat']))
%             save(fullfile(lg_out,img_fn, [img_fn '_cueAlignLick.mat']),  'firstPostRewLickFrame', 'lastPreRewLickFrame', 'tl_rew', 'firstPostRew_lickAlignEvents', 'lastPreRew_lickAlignEvents', 'lickCueAlign', 'lickBurstStart', 'lickCounterVals', 'lickSearch_frames', 'lickDelay_frames', 'lickTC', 'postwin_frames', 'prewin_frames', 'frameRateHz', 'tt', 'ind_early_rew', 'ind_late_rew', 'ind_early_omit', 'ind_late_omit', 'early_rew_time', 'late_rew_time', 'early_omit_time', 'late_omit_time','pct_precue_burst', 'postRew_lickAlignEvents', 'postLick_frames', 'postRew_lickBurstStart','tl', 'ind_prerew_early_rew', 'ind_prerew_late_rew','postRew_lickAlign','preRew_lickBurstHz','postRew_lickBurstHz','ind_low_prerew','ind_high_prerew','ind_low_postrew','ind_high_postrew','ind_low_preomit','ind_high_preomit','ind_low_postomit','ind_high_postomit','HL_lickrate');
%         end
%         if ~exist(fullfile(lg_out,img_fn, [img_fn '_lickResp_preCue.mat']))
%             load(fullfile(jake_dir,img_fn, [img_fn '_lickResp_preCue.mat']))
%             save(fullfile(lg_out,img_fn, [img_fn '_lickResp_preCue.mat']),  'precue_single_lick', 'precue_lick_burst', 'precue_single_lick_df', 'precue_lick_burst_df', 'frameRateHz', 'tl_precue');
%         end
        if ~exist(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            if exist(fullfile(jake_dir,img_fn, [img_fn '_splitImage.mat']))
                load(fullfile(jake_dir,img_fn, [img_fn '_splitImage.mat']))
                save(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']),  'split_img', 'x', 'y','maskCat','maskCat_map');
            end
        end
    end
end