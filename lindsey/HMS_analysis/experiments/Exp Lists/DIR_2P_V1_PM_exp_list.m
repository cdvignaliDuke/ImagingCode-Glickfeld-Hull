exp_list = [];
exp_list.mouse_mat = {'Y13' 'Y13' 'AC44' 'DR9' 'AC42' 'M13' 'M22' 'CM114' 'CM118' 'CM118'};
exp_list.date_mat = {'110509' '110517' '110814' '110922' '111002' '111004' '111201' '120917' '120920' '120921'};
exp_list.runs_mat = {[5] [6] [5:6] [5:6] [1:3] [1:3] [4:6] [1:4] [5:7] [1:4]};
exp_list.prot_mat = {3 3 3 2 1 1 2 1 1 1};
exp_list.run_mat = {1 1 1 1 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1 1 1};
exp_list.SF_mat = {[0.16] [0.16] [0.02 0.08 0.32] [0.16] [0.16 0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04] [0.16 0.04]};
exp_list.TF_mat = {[0.5] [0.5] [1] [1] [2 8] [2 8] [2 8] [2 8] [2 8] [2 8]};
exp_list.dirs_mat = {8 8 8 8 8 8 8 16 16 16};
exp_list.depth_mat = {75 150 75 95 100 80 80 75 75 75};
exp_list.zoom_mat = {1 1 1.5 1.5 1.5 1 1 1 1 1.5};
exp_list.paired_mat = {[] [] [] [1:4] [4:6] [4:5] [1:3] [] [] []};
exp_list.nframes_mat = {12 12 12 12 12 12 12 12 12 12};

fn_out = 'G:\users\lindsey\analysisLG\experiments\DIR_2P_PM\DIR_2P_V1_PM_exp_list.mat';
save(fn_out, 'exp_list');
