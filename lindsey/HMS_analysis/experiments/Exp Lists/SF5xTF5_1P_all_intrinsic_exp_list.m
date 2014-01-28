exp_list = [];
exp_list.mouse_mat = {'M32' 'M38' 'M41' 'M46' 'M47'};
exp_list.date_mat = {'111213' '111214' '120104' '120117' '120120'};
exp_list.runs_mat = {4:7 3:6 1:4 2:5 2:5};
exp_list.prot_mat = {1 1 1 1 1};
exp_list.run_mat = {0 0 0 0 0};
exp_list.blanks_mat = {1 1 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_1P_intrinsic\SF5xTF5_1P_all_intrinsic_exp_list.mat';
save(fn_out, 'exp_list');
