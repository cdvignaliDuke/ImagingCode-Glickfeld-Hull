exp_list = [];
exp_list.mouse_mat = {'VC1' 'VC1' 'VC1' 'CM94' 'CM94' 'CM94' 'CM94'};
exp_list.date_mat = {'120621' '120621' '120623' '120718' '120718' '120723' '120723'};
exp_list.runs_mat = {[1:2] [3:4] [1:2] [1:2] [3:4] [1:3] [4:5]};
exp_list.prot_mat = {1 2 1 1 2 1 2};
exp_list.run_mat = {1 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1};
exp_list.dir_mat = {1 1 1 1 1 1 1};
exp_list.depth_mat = {75 170 75 75 150 75 150};
exp_list.zoom_mat = {1.5 1.5 1.5 1.5 1.5 1.5 1.5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_PM\SF5xTF5_2P_V1L5_PM_exp_list.mat';
save(fn_out, 'exp_list');
