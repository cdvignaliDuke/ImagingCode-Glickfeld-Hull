exp_list = [];
exp_list.mouse_mat = {'AC39' 'AC44' 'Y26' 'AC45' 'Y18' 'Y18' 'DR9' 'M13' 'M14' 'M14' 'M22'};
exp_list.date_mat = {'110809' '110817' '110819' '110820' '110907' '110907' '110915'  '111001' '111123' '111123' '111129'};
exp_list.runs_mat = {[4:6] [1:3] [3:5] [1:4] [5:6] [7:9] [1:3] [4:6] [2:3] [4:5] [1:3]};
exp_list.prot_mat = {1 1 2 1 3 4 1 2 1 2 1};
exp_list.run_mat = {1 1 1 1 1 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1 1 1 1};
exp_list.dir_mat = {2 2 2 2 2 2 2 1 1 1 2};
exp_list.depth_mat = {50 80 60 60 60 110 85 95 80 150 80};
exp_list.zoom_mat = {1 1 1.5 1 1 1 1.5 1.5 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_LM\SF5xTF5_2P_V1_LM_exp_list.mat';
save(fn_out, 'exp_list');
