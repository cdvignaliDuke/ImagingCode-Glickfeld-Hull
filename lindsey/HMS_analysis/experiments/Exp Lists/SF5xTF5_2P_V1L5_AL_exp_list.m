exp_list = [];
exp_list.mouse_mat = {'CM94' 'CM94' 'CM94' 'CM94'};
exp_list.date_mat = {'120804' '120806' '120808' '120808'};
exp_list.runs_mat = {[1:2] [1:2] [1:3] [4:5]};
exp_list.prot_mat = {1 1 1 1};
exp_list.run_mat = {1 1 1 1};
exp_list.blanks_mat = {1 1 1 1};
exp_list.dir_mat = {1 1 1 1};
exp_list.depth_mat = {75 75 75 120};
exp_list.zoom_mat = {1.5 1.5 1.5 1.5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_AL\SF5xTF5_2P_V1L5_AL_exp_list.mat';
save(fn_out, 'exp_list');
