exp_list = [];
exp_list.mouse_mat = {'M11'};
exp_list.date_mat = {'111004'};
exp_list.runs_mat = {[8:11]};
exp_list.prot_mat = {2};
exp_list.run_mat = {1};
exp_list.blanks_mat = {1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_1P_axon\SF5xTF5_1P_AL_axon_exp_list.mat';
save(fn_out, 'exp_list');
