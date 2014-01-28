exp_list = [];
exp_list.mouse_mat = {'Y27' 'Y25' 'Y33' 'M10' 'Y28' 'M21'};
exp_list.date_mat = {'110908'  '110915' '111012' '111116' '111120' '111126'};
exp_list.runs_mat = {2:4 3:7 3:6 3:6 2:5 4:7};
exp_list.prot_mat = {1 1 1 1 1 1};
exp_list.run_mat = {1 1 1 0 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_1P_axon\SF5xTF5_1P_LM_axon_exp_list.mat';
save(fn_out, 'exp_list');
