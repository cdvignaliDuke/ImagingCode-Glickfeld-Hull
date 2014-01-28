exp_list = [];
exp_list.mouse_mat = {'Y25' 'M10' 'Y28' 'M21' 'M23' 'M23'};
exp_list.date_mat = {'111005' '111112' '111125' '111129' '111220' '111220'};
exp_list.runs_mat = {[1:4] [5:7] [1:2] [1:2] [1:2] [3:4]};
exp_list.prot_mat = {1 3 1 1 1 2};
exp_list.run_mat = {1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1};
exp_list.dir_mat = {2 1 1 1 1 1};
exp_list.depth_mat = {70 75 80 90 75 150};
exp_list.zoom_mat = {1 1 1 1 1.5 1.5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_PM\SF5xTF5_2P_LM_PM_exp_list.mat';
save(fn_out, 'exp_list');
