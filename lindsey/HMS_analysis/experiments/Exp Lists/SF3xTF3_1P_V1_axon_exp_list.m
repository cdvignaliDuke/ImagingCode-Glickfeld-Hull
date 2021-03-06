exp_list = [];
exp_list.mouse_mat = {'X32'  'Y13'  'LG29'  'Y2'  'DR8'  'Y14'  'AC39'  'Y26'  'AC45'  'Y18'  'DR9'};
exp_list.date_mat = {'110509'  '110511'  '110515'  '110517'  '110706'  '110706'  '110808'  '110818'  '110818'  '110830'  '110912'};
exp_list.runs_mat = {[2]  [2]  [2]  [4]  [5]  [5 6]  [5 6]  [2 3 4]  [4 5 6]  [4 5 6]  [3 4]};
exp_list.prot_mat = {1 1 1 1 3 2 1 1 1 1 1};
exp_list.run_mat = {0 0 1 1 0 1 1 1 1 1 1};
exp_list.blanks_mat = {0 0 1 1 1 1 1 1 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF3xTF3_1P_axon\SF3xTF3_1P_V1_axon_exp_list.mat';
save(fn_out, 'exp_list');

exp_list = [];
exp_list.mouse_mat = {'X32'  'Y13'  'Y14'  'AC39'  'Y26'  'AC45'  'Y18'  'DR9'  'M13'  'AC42'  'M14'  'M22'};
exp_list.date_mat = {'110509'  '110511'  '110706'  '110808'  '110818'  '110818'  '110830'  '110912' '110927' '110929' '111125' '111126'};
exp_list.runs_mat = {[2]  [2]  [5 6]  [5 6]  [2 3 4]  [4 5 6]  [4 5 6]  [3 4]  [6 7 8] [3 4 5] [3 4 5] [4 5 6]};
exp_list.prot_mat = {1 1 2 1 1 1 1 1 2 1 1 1};
exp_list.run_mat = {0 0 1 1 1 1 1 1 1 1 1 1};
exp_list.blanks_mat = {0 0 1 1 1 1 1 1 1 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF3xTF3_1P_axon\SF3xTF3_1P_V1_axon_exp_list.mat';
save(fn_out, 'exp_list');