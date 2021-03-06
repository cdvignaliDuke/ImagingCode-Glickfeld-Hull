exp_list = [];
exp_list.mouse_mat = {'LG29' 'DR7' 'Y13' 'AC45' 'M13' 'CM114' 'CM118' 'CM118'};
exp_list.date_mat = {'110421' '110425' '110517' '110823' '111002' '120918' '120924' '120924'};
exp_list.runs_mat = {[5] [5:6] [3] [5] [1:3] [1:4] [1:2] [3:5]};
exp_list.prot_mat = {1 1 2 2 1 1 1 2};
exp_list.run_mat = {1 1 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1};
exp_list.SF_mat = {[0.04] [0.02] [0.02] [0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04]};
exp_list.TF_mat = {[8] [8] [8] [8] [2 8] [2 8] [2 8] [2 8]};
exp_list.dirs_mat = {8 8 8 8 8 16 16 16};
exp_list.depth_mat = {[] 125 150 80 100 75 100 50};
exp_list.zoom_mat = {1 1 1 1 1.5 1 1.5 1.5};
exp_list.paired_mat = {[] [] [] [1:4] [4:6] [] [] []};
exp_list.nframes_mat = {15 15 12 12 12 12 12 12};

fn_out = 'G:\users\lindsey\analysisLG\experiments\DIR_2P_AL\DIR_2P_V1_AL_exp_list.mat';
save(fn_out, 'exp_list');