exp_list = [];
exp_list.mouse_mat = {'M31' 'M31' 'M41' 'M41' 'CM114' 'CM114'};
exp_list.date_mat = {'120111' '120111' '120127' '120127' '120920' '120920'};
exp_list.runs_mat = {[1:2] [3:4] [1:2] [3:4] [4:5] [6:7]};
exp_list.prot_mat = {1 2 1 2 1 2};
exp_list.run_mat = {1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1};
exp_list.dir_mat = {1 1 1 1 1 1};
exp_list.depth_mat = {75 150 75 150 75 150};
exp_list.zoom_mat = {1.5 1.5 1 1 1.5 1.5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_AM\SF5xTF5_2P_V1_AM_exp_list.mat';
save(fn_out, 'exp_list');


