exp_list = [];
exp_list.mouse_mat = {'AC39'  'Y26'  'AC45'  'Y18'  'DR9' 'AC42' 'M13' 'M14' 'M22'};
exp_list.date_mat = {'110810'  '110820'  '110819'  '110905'  '110914' '110928' '110928' '111122' '111128'};
exp_list.runs_mat = {[3:6] [2:8] [3:7] [2:4]  [10:12] [2:5] [3:6] [3:6] [3:6]};
exp_list.prot_mat = {1 1 1 1 2 1 1 1 1};
exp_list.run_mat = {1 1 0 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_1P_axon\SF5xTF5_1P_V1_axon_exp_list.mat';
save(fn_out, 'exp_list');
