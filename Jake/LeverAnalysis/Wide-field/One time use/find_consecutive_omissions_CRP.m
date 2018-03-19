%check for # of two omission trials in a row for 2P data
clear
file_info_CRP;
bx_dir = 'Z:\Data\2P_imaging\behavior\';
num_consec_om_cell{1} = [];
for session_num = 1:length(days_post)
    b_data = get_bx_data(bx_dir, days_post{session_num});
    rew_om_bool = cellfun(@isempty, b_data.juiceTimesMsCell); %1=reward omission trials
    rew_om_inx = find(rew_om_bool);
    rew_om_inx
    trial_spacing = diff(rew_om_inx);
    num_consec_om = length(find(trial_spacing==1));
    num_consec_om_cell{session_num} = num_consec_om;
end