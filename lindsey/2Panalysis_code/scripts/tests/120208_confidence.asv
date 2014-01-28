P = 2;
matrix = 'SF5xTF5';
inj = 'V1';
sum_base = 'G:\users\lindsey\analysisLG\experiments';
anal_base = '\\zoloto\bigstorlab\Lindsey\Analysis\120208';
mouse = 'Y13';
date = '110509';
userun = [1:2];

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);

fn_lbub =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_lbub_fits.mat']);
load(fn_lbub);

lbub_diff_TFSF = mean(lbub_diff(:,4:5),2);

figure;
scatter(lbub_diff_TFSF, lbub_fits(:,1,4),3);
xlabel('95% confidence of fits')
ylabel('dF/F')

fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_Y13_confidence_vs_dF.pdf']);
        print(gcf, '-dpdf', fn_out);
        
sum_base = 'G:\users\lindsey\analysisLG\experiments';
fn_summary = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(fn_summary);

fn_sorted =  fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
sorted = readtiff(fn_sorted);

figure;
FOV = mean(sorted,3);
n = all_fits(1).expt(1).n(1);
for iCell = 1:n
    pos = all_fits(1).expt(1).bouton(iCell).pos;
    F = mean(mean(FOV(pos(1)-1:pos(1)+1, pos(2)-1:pos(2)+1),2),1);
    scatter(lbub_diff_TFSF(iCell, :), F , 3, 'b');
    hold on
end
ylim([0 100])
xlabel('95% confidence of fits')
ylabel('F')
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_Y13_confidence_vs_F.pdf']);
        print(gcf, '-dpdf', fn_out);
        
figure;
for iCell = 1:n
    pos = all_fits(1).expt(1).bouton(iCell).pos;
    F = mean(mean(FOV(pos(1)-1:pos(1)+1, pos(2)-1:pos(2)+1),2),1);
    scatter(lbub_fits(iCell,1,4), F , 3, 'b');
    hold on
end
xlabel('dF/F')
ylabel('F')
fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_Y13_dFoverF_vs_F.pdf']);
        print(gcf, '-dpdf', fn_out);
        
[lbub_diff_sort lbub_diff_order] = sort(lbub_diff_TFSF);

plotted = [];
figure;
for iCell = 1:50
    ipix= lbub_diff_order(iCell);
    plotted = [plotted; round(2^lbub_fits(ipix,5,4)*10) round(2^lbub_fits(ipix,4,4)*1000)];
    ind = find(plotted(:,1) == round(2^lbub_fits(ipix,5,4)*10) & plotted(:,2) == round(2^lbub_fits(ipix,4,4)*1000));
    ploterr(2^lbub_fits(ipix,5,4), 2^lbub_fits(ipix,4,4),length(ind)) ;
    hold on
end

figure;
TF_range = [];
SF_range = [];
for iCell = 1:340
    ipix= lbub_diff_order(iCell);
    TF_range = [TF_range; 2^lbub_fits(ipix,5,4)-2^lbub_fits(ipix,5,1) 2^lbub_fits(ipix,5,4)-2^lbub_fits(ipix,5,2)];
    SF_range = [SF_range; 2^lbub_fits(ipix,4,4)-2^lbub_fits(ipix,4,1) 2^lbub_fits(ipix,4,4)-2^lbub_fits(ipix,4,2)];
    for iTF = 1:size(TF_range,1)
        
            
    ploterr(2^lbub_fits(ipix,5,4), 2^lbub_fits(ipix,4,4),{2^lbub_fits(ipix,5,1), 2^lbub_fits(ipix,5,2)},{2^lbub_fits(ipix,4,1),2^lbub_fits(ipix,4,2)}) ;
    hold on
end
ind = find(lbub_diff_sort>2);


