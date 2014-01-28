P = 2;
inj = 'V1';
matrix = 'SF5xTF5';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);

areas = ['PM'; 'LM'; 'AL'; 'RL'; 'AM'];

for iArea = 1:5
    TF = [];
    SF = [];
    sigma_TF = [];
    sigma_SF = [];
    dF = [];
    speed = [];
    xi = [];
    plotfit = [];
    data = [];
    nexp = all_fits(iArea).nexp;
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    for iexp = 1:nexp
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        base = 'G:\users\lindsey\analysisLG\active mice';    
        outDir = fullfile(base, mouse,date);
        n = all_fits(iArea).expt(iexp).n(1);
        fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
        load(fn_out);
        for iCell = 1:n
            if lbub_diff(iCell,4)<2;
                if lbub_diff(iCell,5)<2;
                    TF = [TF; all_fits(iArea).expt(iexp).bouton(iCell).TF_fit];
                    SF = [SF; all_fits(iArea).expt(iexp).bouton(iCell).SF_fit];
                    sigma_TF = [sigma_TF; all_fits(iArea).expt(iexp).bouton(iCell).sigma_TF];
                    sigma_SF = [sigma_SF; all_fits(iArea).expt(iexp).bouton(iCell).sigma_SF];
                    dF = [dF; all_fits(iArea).expt(iexp).bouton(iCell).dF_fit];
                    speed = [speed; all_fits(iArea).expt(iexp).bouton(iCell).speed];
                    xi = [xi; all_fits(iArea).expt(iexp).bouton(iCell).xi_fit];
                    plotfit = [plotfit; reshape(all_fits(iArea).expt(iexp).bouton(iCell).plotfit', 1, 25)];
                    data = [data; reshape(all_fits(iArea).expt(iexp).bouton(iCell).dFoverF', 1, 25)];
                end
            end
        end
    end
    Goodfits(iArea).TF = TF;
    Goodfits(iArea).SF = SF;
    Goodfits(iArea).sigma_TF = sigma_TF;
    Goodfits(iArea).sigma_SF = sigma_SF;
    Goodfits(iArea).dF = dF;
    Goodfits(iArea).speed = speed;
    Goodfits(iArea).xi = xi;
    Goodfits(iArea).plotfit = plotfit;
    Goodfits(iArea).data = data;
    Goodfits(iArea).expt = iexp;
end
fn_good = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'good_fits.mat');
save(fn_good, 'Goodfits');