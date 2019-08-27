% find movement onset latency in response to the cue/reward
% used in piezo_analysis after piezo_mov_onset
function trial_onset_lats = piezo_onset_lat_in_trial(mov_inds, reward_omit_inx, unexp_rew_inx, rew_inds, cue_inds);

%set criteria for looking for onset latency.
%no suprathreshold movement 50ms before to 150ms after cue onset.


%group variables together
trial_onset_lats.trial_onset_lats_NR; 
trial_onset_lats.trial_onset_lats_OR;
trial_onset_lats.trial_onset_lats_UR;

return