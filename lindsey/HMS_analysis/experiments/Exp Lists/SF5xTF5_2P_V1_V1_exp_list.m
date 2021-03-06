exp_list = [];
exp_list.mouse_mat = {'CM63' 'CM63' 'CM63' 'CM63' 'VC1' 'VC1' 'VC1' 'CM69' 'CM69' 'CM69' 'CM94' 'CM94' 'CM94'};
exp_list.date_mat = {'120524' '120525' '120527' '120529' '120602' '120603' '120604' '120603' '120604' '120605' '120705' '120709' '120710'};
exp_list.runs_mat = {[1:4] [1:4] [1:4] [1:4] [1:3] [1:3] [4:6] [1:3] [1:3] [1:3] [1:3] [1:3] [1:3]};
exp_list.prot_mat = {1 1 1 1 1 1 2 1 1 1 1 1 1};
exp_list.run_mat = {1 1 1 1 1 1 1 1 1 1 1 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1 1 1 1 1 1};
exp_list.dir_mat = {1 1 1 1 1 1 1 1 1 1 1 1 1};
exp_list.depth_mat = {350 425 150 250 215 290 375 300 215 120 375 320 450};
exp_list.zoom_mat = {5 5 3 4 5 5 5 5 5 2.5 5 5 5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_V1\SF5xTF5_2P_V1_V1_exp_list.mat';
save(fn_out, 'exp_list');
