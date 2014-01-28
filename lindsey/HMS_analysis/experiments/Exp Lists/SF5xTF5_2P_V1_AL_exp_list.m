exp_list = [];
exp_list.mouse_mat = {'DR7' 'Y13' 'X32' 'AC39' 'AC44' 'AC45' 'AC45' 'Y26' 'DR9' 'M13' 'M14' 'M22' 'M31' 'M31'};
exp_list.date_mat = {'110505' '110509' '110512' '110812' '110817' '110823' '110823' '110823' '110920'  '111002' '111126' '111203' '120103' '120103'};
exp_list.runs_mat = {[1] [3:4] [1:2] [1:4] [4:6] [1:4] [6:9] [3:6] [1:4] [4:6] [1:4] [1:4] [1:2] [3:4]};
exp_list.prot_mat = {1 2 1 1 2 1 3 1 1 2 1 1 1 2};
exp_list.run_mat = {0 0 0 1 1 1 1 1 1 1 1 0 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1 1 1 1 1 1 1};
exp_list.dir_mat = {1 1 1 2 2 2 2 2 2 1 2 2 1 1};
exp_list.depth_mat = {70 75 85 90 80 80 50 70 90 100 80 90 75 150};
exp_list.zoom_mat = {1 1 1 1 1 1 1 1.5 1 1.5 1 1 1.5 1.5};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_2P_AL\SF5xTF5_2P_V1_AL_exp_list.mat';
save(fn_out, 'exp_list');


