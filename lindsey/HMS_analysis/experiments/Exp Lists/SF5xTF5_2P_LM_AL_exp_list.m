exp_list = [];
exp_list.mouse_mat = {'M10' 'Y28' 'M23' 'M23' 'M15' 'M15'};
exp_list.date_mat = {'111112' '111126' '111216' '111216' '111219' '111219'};
exp_list.runs_mat = {[1:2] [1:2] [1:2] [3] [1:2] [3:4]};
exp_list.prot_mat = {1 1 1 2 1 2};
exp_list.run_mat = {1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1};
exp_list.dir_mat = {1 1 1 1 1 1};
exp_list.depth_mat = {95 75 75 175 75 150};
exp_list.zoom_mat = {1 1.5 1.5 1.5 1.5 1.5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_AL\SF5xTF5_2P_LM_AL_exp_list.mat';
save(fn_out, 'exp_list');
