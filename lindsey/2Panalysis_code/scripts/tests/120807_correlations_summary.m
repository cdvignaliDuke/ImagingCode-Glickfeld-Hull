P = 2;
inj = 'V1';
matrix = 'SF5xTF5';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);

areas = ['PM'; 'LM'; 'AL'];
Corr = [];
for iArea = 1:3
    nexp = all_fits(iArea).nexp;
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    real = [];
    shuf = [];
    sub = [];
    blank = [];
    for iexp = 1:nexp
        if all_fits(iArea).expt(iexp).n(2)>0
            mouse = char(exp_list.mouse_mat{iexp});
            date = char(exp_list.date_mat{iexp});
            userun = exp_list.runs_mat{iexp};
            base = 'G:\users\lindsey\analysisLG\active mice';    
            outDir = fullfile(base, mouse,date);
            fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_bouton_corr.mat']);
            load(fn_out)
            real = [real; r_real_avg];
            shuf = [shuf; r_shuf_avg];
            sub = [sub; r_sub_avg];
            fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_bouton_corr_blank.mat']);
            load(fn_out)
            blank = [blank; r_blank_avg];
        end
    end
    Corr(iArea).real = real;
    Corr(iArea).sub = sub;
    Corr(iArea).shuf = shuf;
    Corr(iArea).blank = blank;
end