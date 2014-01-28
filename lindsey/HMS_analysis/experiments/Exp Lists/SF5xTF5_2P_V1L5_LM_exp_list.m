exp_list = [];
exp_list.mouse_mat = {'VC1' 'VC1' 'CM94' 'CM94'};
exp_list.date_mat = {'120622' '120622' '120720' '120720'};
exp_list.runs_mat = {[1:2] [3:4] [1:2] [3:4]};
exp_list.prot_mat = {1 2 1 2};
exp_list.run_mat = {1 1 1 1};
exp_list.blanks_mat = {1 1 1 1};
exp_list.dir_mat = {1 1 1 1};
exp_list.depth_mat = {60 150 75 150};
exp_list.zoom_mat = {1.5 1.5 1.5 1.5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_LM\SF5xTF5_2P_V1L5_LM_exp_list.mat';
save(fn_out, 'exp_list');
