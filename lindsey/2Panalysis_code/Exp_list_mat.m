exp_list = [];
exp_list.mouse_mat = {'AC44' 'AC45' 'Y26' 'DR9' 'M13' 'AC42'};
exp_list.date_mat = {'110802' '110804' '110818' '110913' '110924' '110925'};
exp_list.runs_mat = {2:4 3:5 5:8 2:4 4:6 3:6};
exp_list.prot_mat = {1 1 2 1 2 1};
exp_list.run_mat = {1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_1P\SF5xTF5_V1_exp_list.mat';
save(fn_out, 'exp_list');