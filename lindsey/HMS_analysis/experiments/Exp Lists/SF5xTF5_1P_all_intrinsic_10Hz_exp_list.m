exp_list = [];
exp_list.mouse_mat = {'M38' 'M32' 'M41' 'M46' 'M47'};
exp_list.date_mat = {'111215' '111218' '120103' '120116' '120126'};
exp_list.runs_mat = {3:6 1:5 1:2 2:3 2:3};
exp_list.prot_mat = {2 1 1 1 1};
exp_list.run_mat = {0 0 0 0 0};
exp_list.blanks_mat = {1 1 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_1P_intrinsic_10Hz\SF5xTF5_1P_all_intrinsic_10Hz_exp_list.mat';
save(fn_out, 'exp_list');
