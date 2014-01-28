P = 1;
matrix = 'SF5xTF5';
inj = 'all';
image = 'intrinsic';
nCond = 25;
nON=10;
nOFF=10;
nPlanes = 1;
begin = 7;
TFSFetc = [1:2];
pre_win = [1 4];
post_win = [5 14];
sum_base = 'G:\users\lindsey\analysisLG\experiments';

list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);

for iexp = 1:nexp
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_protocol = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    
    base = 'G:\users\lindsey\analysisLG\active mice';
    running_base = 'G:\users\lindsey\dataLG\Running data';
    outDir = fullfile(base, mouse,date);
    
    if run==1
    Running_LG
    end
end
     %load roi_avg
    fn_avg = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_roi_avg.mat']);
    load(fn_avg);

    load(reps)
    load
    roi_sort_norm
end