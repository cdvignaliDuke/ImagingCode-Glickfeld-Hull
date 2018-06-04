% verify_trial_types_with_RTs
%use the indeces of the trials used for imaging to verify that tha trial
%type classifications are correct. Should be run in main2


load([dest 'parse_behavior.mat']);
b_data = input;

%correct trials
corr_trials_bdata_inx = find(trial_outcome.corr_inx);
corr_trials_bdata_inx([trial_outcome.succ_rm_event_idx]) = [];
RT_in_ms = cell2mat( b_data.reactTimesMs(corr_trials_bdata_inx));
assert(min(RT_in_ms)>=200);

%early trials
fail_trials_bdata_inx = find(trial_outcome.early_inx);
fail_trials_bdata_inx([trial_outcome.fail_rm_event_idx]) = [];
RT_in_ms = cell2mat(b_data.reactTimesMs([fail_trials_bdata_inx]));
assert(max(RT_in_ms)<1);

%too fast corrects
tf_trials_bdata_inx = find(trial_outcome.tooFast_inx);
if ~isempty(find(trial_outcome.tooFast_inx))
    tf_trials_bdata_inx([trial_outcome.TF_rm_event_idx]) = [];
    RT_in_ms = cell2mat(b_data.reactTimesMs([tf_trials_bdata_inx]));
    assert(max(RT_in_ms)<200 & min(RT_in_ms)>0);
end
