exp_list = [];
exp_list.mouse_mat = {'M15' 'M23' 'M23'};
exp_list.date_mat = {'111221' '111221' '111221'};
exp_list.runs_mat = {[1:2] [3] [4:5]};
exp_list.prot_mat = {1 2 3};
exp_list.run_mat = {1 1 1};
exp_list.blanks_mat = {1 1 1};
exp_list.dir_mat = {1 1 1};
exp_list.depth_mat = {75 75 95};


fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_AM\SF5xTF5_2P_LM_AM_exp_list.mat';
save(fn_out, 'exp_list');
