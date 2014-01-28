clear all
areas = ['PM'; 'LM'; 'AL'];
inj = 'V1';
P = 2;
nPlanes = 1;
begin = 1;
Nshuf = 500;

for iArea = 1:3
    matrix = 'DIR';
    image = areas(iArea,:);
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
        SFs = exp_list.SF_mat{iexp};
        TFs = exp_list.TF_mat{iexp};
        dirs = exp_list.dirs_mat{iexp};
        pair_run = exp_list.paired_mat{iexp};
        nframes = exp_list.nframes_mat{iexp};
        
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        fn_local= fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_local_max.mat']);
        load(fn_local);
        
        figure;
        imagesq(local_max);
        title([mouse ' ' date]) 
    end
end