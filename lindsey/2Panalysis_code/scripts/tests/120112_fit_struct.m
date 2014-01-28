P = 2;
matrix = 'SF5xTF5';
image = 'PM';
inj = 'V1';
nON=12;
nOFF=12;
Nshuf = 500;
SF_vec0 = [.32 .16 .08 .04 .02]; %flipped to have low to high SF in square  %flipud
TF_vec0 = [1 2 4 8 15];

sum_base = 'G:\users\lindsey\analysisLG\experiments';
base = 'G:\users\lindsey\analysisLG\active mice';

list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp = size(exp_list.mouse_mat,2);
mouse_list = [];
mice = exp_list.mouse_mat;

iexp = 3;
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_prot = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    dir = exp_list.dir_mat{iexp};
    
    base = 'G:\users\lindsey\analysisLG\active mice';    
    outDir = fullfile(base, mouse,date);

fn_local= fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_local_max.mat']);
load(fn_local);
    fn_fit = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_Fit_struct.mat']);   
load(fn_fit);
fn_lbub = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
load(fn_lbub);

figure;
start = 1;
fig = 1;
for iCell = 1:n_pix    
    if start>64
        fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120112', [date '_' mouse '_run' num2str(userun) '_Data_vs_Fit_' num2str(fig) '.pdf']);
        print(gcf, '-dpdf', fn_out);
        figure;
        fig = 1+fig;
        start = 1;
    end
    subplot(8, 8, start);
    imagesq(Fit_struct(iCell).True.s_.data);
    subplot(8, 8, start+1);
    imagesq(Fit_struct(iCell).True.s_.k2b_plot);
    if lbub_diff(iCell,4)<2 
        if lbub_diff(iCell,5)<2
            title('**');
        end
    end
    colormap(gray);
    start= start+2;
    if iCell == n_pix
        fn_out = fullfile('\\Zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging\Analysis\120112', [date '_' mouse '_run' num2str(userun) '_Data_vs_Fit_' num2str(fig) '.pdf']);
        print(gcf, '-dpdf', fn_out);
    end
end