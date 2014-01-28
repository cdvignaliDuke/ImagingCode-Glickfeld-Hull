areas = ['PM'; 'LM'; 'AL'];
CM_mat_avg = zeros(3,2,3);
CM_mat_sem = zeros(3,2,3);
resp_avg_avg = zeros(3,25,3);
resp_norm_avg = zeros(3,25,3);
nexp = zeros(3,1);
for iArea = 1:3;
P = 2;
matrix = 'SF5xTF5';
image = areas(iArea,:);
inj = 'V1';
nCond = 25;
nON=12;
nOFF=12;

sum_base = 'G:\users\lindsey\analysisLG\experiments';
base = 'G:\users\lindsey\analysisLG\active mice';

list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
load(list_fn);
nexp(iArea,:) = size(exp_list.mouse_mat,2);
mouse_list = [];
mice = exp_list.mouse_mat;

CM_mat_all = zeros(3, 2, nexp);
resp_avg_all = zeros(3,25,nexp);
resp_norm_all = zeros(3,25,nexp);

for iexp = 1:nexp
    mouse = char(exp_list.mouse_mat{iexp});
    date = char(exp_list.date_mat{iexp});
    userun = exp_list.runs_mat{iexp};
    count_prot = exp_list.prot_mat{iexp};
    run = exp_list.run_mat{iexp};
    blanks = exp_list.blanks_mat{iexp};
    
    outDir = fullfile(base, mouse,date);
    
    fn_CM = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_CM_data.mat']);
    load(fn_CM);
    fn_resp = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp_avg.mat']);
    load(fn_resp);
    
    CM_mat_all(:,:,iexp) = CM_data;
    resp_avg_all(:,:,iexp) = resp_avg;
    for ithresh = 1:3
        resp_norm_all(ithresh,:,iexp) = resp_avg(ithresh,:)./squeeze(max(resp_avg(ithresh,:),[],2));
    end
end

CM_mat_avg(:,:,iArea) = nanmean(CM_mat_all,3);
CM_mat_sem(:,:,iArea) = nanstd(CM_mat_all,[],3)/sqrt(nexp);

resp_avg_avg(:,:,iArea) = nanmean(resp_avg_all,3);
resp_norm_avg(:,:,iArea) = nanmean(resp_norm_all,3);
end

figure;
col = ['k' 'b' 'r'];
for ithresh = 1:3
    h= subplot(2,2,1);
    errorbar(1:3, CM_mat_avg(ithresh,1,:), CM_mat_sem(ithresh,1,:), col(ithresh));
    ylim([0 .1])
    ylabel('Spatial Frequency (cpd)')
    set(h,'XTick', [1:length(areas)]);    
    set(h,'XTickLabel', areas);
    hold on
    h= subplot(2,2,2);
    errorbar(1:3, CM_mat_avg(ithresh,2,:), CM_mat_sem(ithresh,2,:), col(ithresh));
    ylim([0 7])
    ylabel('Temporal Frequency (Hz)')
    set(h,'XTick', [1:length(areas)]);
    set(h,'XTickLabel', areas);
    hold on
end

figure;
start = 1;
for ithresh = 1:3
    for iArea= 1:3
        subplot(3,4,start)
        imagesq(reshape(resp_avg_avg(ithresh,:,iArea), [5 5])');
        colormap gray
        start = start+1;
    end
end
subplot(3,3,1)
ylabel('1%')
subplot(3,3,4)
ylabel('10%')
subplot(3,3,7)
ylabel('ALL')

subplot(3,3,1)
title('PM')
subplot(3,3,2)
title('LM')
subplot(3,3,3)
title('AL')

suptitle('Average')
suptitle('Normalized')



figure;
start = 1;
for ithresh = 1:3
    for iArea= 1:3
        subplot(3,3,start)
        imagesq(reshape(resp_norm_avg(ithresh,:,iArea), [5 5])');
        colormap gray
        start = start+1;
    end
end
        