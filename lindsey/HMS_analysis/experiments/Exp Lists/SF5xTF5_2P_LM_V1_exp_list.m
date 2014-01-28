exp_list = [];
exp_list.mouse_mat = {'Y27' 'Y27' 'M10' 'M10' 'Y28' 'M21' 'M15'};
exp_list.date_mat = {'110914' '110916' '111111' '111118' '111119' '111203' '111222'};
exp_list.runs_mat = {[3] [3:5] [1:2] [1:4] [1:2] [1:2] [1:2]};
exp_list.prot_mat = {2 2 1 1 1 1 1};
exp_list.run_mat = {1 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1};
exp_list.dir_mat = {2 2 1 2 1 1 1};
exp_list.depth_mat = {70 75 80 80 80 60 75};
exp_list.zoom_mat = {1.5 1.5 1.5 1 1 1 1.5};


fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_V1\SF5xTF5_2P_LM_V1_exp_list.mat';
save(fn_out, 'exp_list');

