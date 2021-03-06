exp_list = [];
exp_list.mouse_mat = {'Y11' 'DR4' 'AC44' 'AC45' 'Y26' 'DR9' 'M13' 'AC42' 'Y27' 'Y25' 'M11' 'Y32' 'Y33' 'M17' 'DR10' 'M10' 'Y28' 'M22' 'M14' 'M21' 'M15' 'M23'};
exp_list.date_mat = {'110525' '110709' '110802' '110804' '110818' '110913' '110924' '110925' '110905'  '110912' '111001' '111010' '111010' '111024' '111025' '111026' '111108' '111111' '111112' '111118' '111201' '111202'};
exp_list.runs_mat = {1:4 4:6 2:4 3:5 5:8 2:4 4:6 3:6 3:6 3:5 5:8 2:5 2:5 2:5 5:8 3:6 3:6 5:8 5:7 4:7 2:4 4:6};
exp_list.prot_mat = {1 1 1 1 2 1 2 1 1 1 1 1 1 1 1 1 1 2 2 1 1 1};
exp_list.run_mat = {1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 1 1};
exp_list.blanks_mat = {1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1};

fn_out = 'G:\users\lindsey\analysisLG\experiments\SF5xTF5_1P_soma\SF5xTF5_1P_all_soma_exp_list.mat';
save(fn_out, 'exp_list');