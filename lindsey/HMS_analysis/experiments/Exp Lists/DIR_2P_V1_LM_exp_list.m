exp_list = [];
exp_list.mouse_mat = {'Y26' 'Y26' 'Y18' 'Y18' 'DR9' 'M13' 'AC42' 'M22' 'CM114' 'CM118'};
exp_list.date_mat = {'110824' '110824' '110907' '110907' '110915' '111001' '111001' '111129' '120914' '120919'};
exp_list.runs_mat = {[1:3] [4:6] [1:2] [3:4] [4:5] [1:3] [1:2] [4:5] [1:4] [1:3]};
exp_list.prot_mat = {1 2 1 2 2 1 1 2 1 1};
exp_list.run_mat = {1 1 1 1 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1 1 1};
exp_list.SF_mat = {[0.04] [0.04] [0.04] [0.04] [0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04]};
exp_list.TF_mat = {[2] [2] [2] [2] [2] [2 8] [2 8] [2 8] [2 8] [2 8]};
exp_list.dirs_mat = {8 8 8 8 8 8 8 8 16 16};
exp_list.depth_mat = {50 120 50 110 85 95 90 80 75 75};
exp_list.zoom_mat = {1.5 1.5 1 1 1.5 1.5 1 1 1 1.5};
exp_list.paired_mat = {[] [] [] [5:6] [1:3] [4:6] [] [1:3] [] []};
exp_list.nframes_mat = {12 12 12 12 12 12 12 12 12 12};

fn_out = 'G:\users\lindsey\analysisLG\experiments\DIR_2P_LM\DIR_2P_V1_LM_exp_list.mat';
save(fn_out, 'exp_list');
