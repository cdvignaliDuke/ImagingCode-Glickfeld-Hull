close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

sigCells = [];
fractSigCells = [];
respCells = [];
resp_ind_all = [];

p_anova_all_all = [];
b_all_all = [];
amp_all_all = [];
pha_all_all = [];
per_all_all = [];
sse_all_all = [];
R_square_all_all = [];
yfit_all_all = [];

p_anova_shuf_all = [];
b_shuf_all = [];
amp_shuf_all = [];
pha_shuf_all = [];
per_shuf_all = [];
sse_shuf_all = [];
R_square_shuf_all = [];
yfit_shuf_all = [];

p_anova_downsamp_all = [];
b_downsamp_all = [];
amp_downsamp_all = [];
pha_downsamp_all = [];
per_downsamp_all = [];
sse_downsamp_all = [];
R_square_downsamp_all = [];
yfit_downsamp_all = [];

plaidSI_all = [];
testPI_all = [];
OSI_all = [];
max_dir_all = [];

ds = ['CrossOriRandPhase_ExptList'];
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);

totCells = zeros(nexp,1);
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).coFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);

    fprintf([mouse ' ' date '\n'])

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_phaseFits.mat']))

    nCells = size(resp_cell{end,end,end},1);
    totCells(iexp,:) = nCells;

    p_anova_all(find(p_anova_all==0)) = NaN;
    p_anova_all_all = [p_anova_all_all;  p_anova_all];
    b_all_all = [b_all_all; b_hat_all];
    amp_all_all = [amp_all_all; amp_hat_all];
    pha_all_all = [pha_all_all; pha_hat_all];
    per_all_all = [per_all_all; per_hat_all];
    R_square_all_all = [R_square_all_all; R_square_all];
    sse_all_all = [sse_all_all; sse_all];

    p_anova_shuf(find(p_anova_shuf==0)) = NaN;
    p_anova_shuf_all = [p_anova_shuf_all;  p_anova_shuf];
    b_shuf_all = [b_shuf_all; b_hat_shuf];
    amp_shuf_all = [amp_shuf_all; amp_hat_shuf];
    pha_shuf_all = [pha_shuf_all; pha_hat_shuf];
    per_shuf_all = [per_shuf_all; per_hat_shuf];
    R_square_shuf_all = [R_square_shuf_all; R_square_shuf];
    sse_shuf_all = [sse_shuf_all; sse_shuf];
    
    p_anova_downsamp(find(p_anova_downsamp==0)) = NaN;
    p_anova_downsamp_all = [p_anova_downsamp_all;  p_anova_downsamp];
    b_downsamp_all = [b_downsamp_all; b_hat_downsamp];
    amp_downsamp_all = [amp_downsamp_all; amp_hat_downsamp];
    pha_downsamp_all = [pha_downsamp_all; pha_hat_downsamp];
    per_downsamp_all = [per_downsamp_all; per_hat_downsamp];
    R_square_downsamp_all = [R_square_downsamp_all; R_square_downsamp];
    sse_downsamp_all = [sse_downsamp_all; sse_downsamp];

    yfit_all_all = cat(1,yfit_all_all, yfit_all);
    yfit_shuf_all = cat(1,yfit_shuf_all, yfit_shuf);
    yfit_downsamp_all = cat(1,yfit_downsamp_all, yfit_downsamp);

    resp_ind_all = [resp_ind_all; resp_ind+sum(totCells(1:iexp-1,:),1)];

    plaid_resp = mean(resp_cell{end,end,1},2);
    mask_resp = mean(resp_cell{end,1,1},2);
    test_resp = mean(resp_cell{1,end,1},2);
    plaid_resp(find(plaid_resp<0)) = 0;
    mask_resp(find(mask_resp<0)) = 0;
    test_resp(find(test_resp<0)) = 0;

    plaidSI_all = [plaidSI_all; (plaid_resp-(mask_resp+test_resp)) ./ (plaid_resp + mask_resp + test_resp)];
    testPI_all = [testPI_all; abs((test_resp-mask_resp) ./ (mask_resp+test_resp))];
    
    ImgFolder = expt(iexp).dirFolder;
    nrun = length(ImgFolder);
    run_str = catRunName(cell2mat(ImgFolder), nrun);
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    OSI_all = [OSI_all; OSI_resp];
    max_dir_all = [max_dir_all; max_ind_dir];
end

%%
figure;
subplot(2,2,1)
scatter(p_anova_all_all(resp_ind_all,:), amp_all_all(resp_ind_all,:))
hold on
xlabel('anova p-value')
ylabel('Sine amplitude')
ylim([0 2])

figure;
subplot(2,2,1)
Rsq_temp = R_square_all_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all,:))
n = sum(~isnan(Rsq_temp(resp_ind_all,:)));
hold on
Rsq_temp = R_square_shuf_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all,:))
Rsq_temp = R_square_downsamp_all;
Rsq_temp(find(Rsq_temp<0)) = 0;
cdfplot(Rsq_temp(resp_ind_all,:))
xlabel('Rsquared')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled','Downsampled'},'location','southeast')
subplot(2,2,2)
cdfplot(amp_all_all(resp_ind_all,:))
hold on
cdfplot(amp_shuf_all(resp_ind_all,:))
cdfplot(amp_downsamp_all(resp_ind_all,:))
xlabel('Sine Amplitude')
xlim([0 1])
legend({['Data- n = ' num2str(n)],'Shuffled','Downsampled'},'location','southeast')
print(fullfile(summaryDir, 'randPhase_RsqAmp_AllTrialFits_8phase.pdf'),'-dpdf','-fillpage')


figure;
subplot(2,2,1)
ind = resp_ind_all;
cdfplot(amp_all_all(ind(find(testPI_all(ind,:)<0.3),:),1))
hold on
cdfplot(amp_all_all(ind(find(testPI_all(ind,:)>0.9),:),1))
xlabel('Sine amplitude')
title('')
legend({['PI<0.3- n=' num2str(length(find(testPI_all(ind,:)<0.3)))] ,['PI>0.9- n=' num2str(length(find(testPI_all(ind,:)>0.9)))]},'location','southeast')
subplot(2,2,3)
cdfplot(plaidSI_all(ind(find(testPI_all(ind,:)<0.3),:),1))
hold on
cdfplot(plaidSI_all(ind(find(testPI_all(ind,:)>0.9),:),1))
xlabel('Suppression index')
title('')
legend({'PI<0.3','PI>0.9'},'location','northwest')

subplot(2,2,2)
ind = resp_ind_all;
cdfplot(R_square_all_all(ind(find(testPI_all(ind,:)<0.3),:),1))
hold on
cdfplot(R_square_all_all(ind(find(testPI_all(ind,:)>0.9),:),1))
xlabel('Rsq')
title('')
legend({'PI<0.3', 'PI>0.9'},'location','southeast')

subplot(2,2,4)
ind = resp_ind_all;
cdfplot(amp_all_all(ind(find(plaidSI_all(ind,:)<0),:),1))
hold on
cdfplot(amp_all_all(ind(find(plaidSI_all(ind,:)>0),:),1))
xlabel('Sine Amp')
title('')
legend({['SI<0- n=' num2str(length(find(plaidSI_all(ind,:)<0)))] ,['SI>0- n=' num2str(length(find(plaidSI_all(ind,:)>0)))]},'location','southeast')
print(fullfile(summaryDir, 'randPhase_sineAmpvsPrefIndex_8phase.pdf'),'-dpdf','-fillpage')

figure;
subplot(2,2,1)
ind = resp_ind_all;
cdfplot(amp_all_all(ind(find(OSI_all(ind,:)<0.5),:),1))
hold on
cdfplot(amp_all_all(ind(find(OSI_all(ind,:)>0.5),:),1))
xlabel('Sine amplitude')
title('')
legend({['OSI<0.5- n=' num2str(length(find(OSI_all(ind,:)<0.5)))] ,['OSI>0.5- n=' num2str(length(find(OSI_all(ind,:)>0.5)))]},'location','southeast')
subplot(2,2,3)
cdfplot(plaidSI_all(ind(find(OSI_all(ind,:)<0.5),:),1))
hold on
cdfplot(plaidSI_all(ind(find(OSI_all(ind,:)>0.5),:),1))
xlabel('Suppression index')
title('')
legend({'OSI<0.5','OSI>0.5'},'location','northwest')

subplot(2,2,2)
ind = resp_ind_all;
cdfplot(R_square_all_all(ind(find(OSI_all(ind,:)<0.5),:),1))
hold on
cdfplot(R_square_all_all(ind(find(OSI_all(ind,:)>0.5),:),1))
xlabel('Rsq')
title('')
legend({'OSI<0.5', 'OSI>0.5'},'location','southeast')
print(fullfile(summaryDir, 'randPhase_sineAmpvsOSI_8phase.pdf'),'-dpdf','-fillpage')

figure;
[n edges bin] = histcounts(testPI_all);
for i = 1:length(n)
    ind = intersect(find(R_square_all_all>0.1),intersect(resp_ind_all, find(bin == i)));
    subplot(2,2,1)
    errorbar(nanmean(testPI_all(ind,:),1), nanmean(amp_all_all(ind,:),1),nanstd(amp_all_all(ind,:),[],1)./sqrt(length(ind)),nanstd(amp_all_all(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
    subplot(2,2,2)
    errorbar(nanmean(testPI_all(ind,:),1), nanmean(plaidSI_all(ind,:),1),nanstd(plaidSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(plaidSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(testPI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
end
subplot(2,2,1)
ylim([0 1])
ylabel('Sine Amplitude')
xlabel('Ori preference index')
subplot(2,2,2)
ylim([-0.5 .5])
ylabel('Suppression index')
xlabel('Ori preference index')
[n edges bin] = histcounts(OSI_all,5);
for i = 1:length(n)
    ind = intersect(find(R_square_all_all>0.1),intersect(resp_ind_all, find(bin == i)));
    subplot(2,2,3)
    errorbar(nanmean(OSI_all(ind,:),1), nanmean(amp_all_all(ind,:),1),nanstd(amp_all_all(ind,:),[],1)./sqrt(length(ind)),nanstd(amp_all_all(ind,:),[],1)./sqrt(length(ind)),nanstd(OSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(OSI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
    subplot(2,2,4)
    errorbar(nanmean(OSI_all(ind,:),1), nanmean(plaidSI_all(ind,:),1),nanstd(plaidSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(plaidSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(OSI_all(ind,:),[],1)./sqrt(length(ind)),nanstd(OSI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
end
subplot(2,2,3)
ylim([0 1])
ylabel('Sine Amplitude')
xlabel('OSI')
subplot(2,2,4)
ylim([-0.5 .5])
ylabel('Suppression index')
xlabel('OSI')
print(fullfile(summaryDir, 'randPhase_OSIorPrefIndexBins_8phase.pdf'),'-dpdf','-fillpage')


figure;
ind = resp_ind_all;
max_dir_ind_reset = max_dir_all;
max_dir_ind_reset(find(max_dir_all>180)) = max_dir_ind_reset(find(max_dir_all>180))-360;
pref_dir_ind = intersect(ind,intersect(find(max_dir_ind_reset>-22.5),find(max_dir_ind_reset<22.5)));
subplot(2,2,1)
scatter(testPI_all(ind,:), amp_all_all(ind,:),'ok')
hold on
scatter(testPI_all(pref_dir_ind,:), amp_all_all(pref_dir_ind,:),'or')
ylim([0 1])
xlim([0 1])
ylabel('Sine Amplitude')
xlabel('Stim preference index')
% subplot(2,3,2)
% scatter(testPI_all(ind,:), amp_all_all(ind,:)-amp_shuf_all(ind,:),'ok')
% hold on
% scatter(testPI_all(pref_dir_ind,:), amp_all_all(pref_dir_ind,:)-amp_shuf_all(pref_dir_ind,:),'or')
% ylim([0 1])
% xlim([0 1])
% ylabel('Sine Amplitude (Data-Shuf)')
% xlabel('Stim preference index')
subplot(2,2,2)
scatter(testPI_all(ind,:), plaidSI_all(ind,:),'ok')
hold on
scatter(testPI_all(pref_dir_ind,:), plaidSI_all(pref_dir_ind,:),'or')
ylim([-1 1])
xlim([0 1])
ylabel('Suppression index')
xlabel('Stim preference index')
subplot(2,2,3)
scatter(OSI_all(ind,:), amp_all_all(ind,:),'ok')
hold on
scatter(OSI_all(pref_dir_ind,:), amp_all_all(pref_dir_ind,:),'or')
ylim([0 1])
xlim([0 1])
ylabel('Sine Amplitude')
xlabel('OSI')
% subplot(2,3,5)
% scatter(OSI_all(ind,:), amp_all_all(ind,:)-amp_shuf_all(ind,:),'ok')
% hold on
% scatter(OSI_all(pref_dir_ind,:), amp_all_all(pref_dir_ind,:)-amp_shuf_all(pref_dir_ind,:),'or')
% ylim([0 1])
% xlim([0 1])
% ylabel('Sine Amplitude')
% xlabel('OSI')
subplot(2,2,4)
scatter(OSI_all(ind,:), plaidSI_all(ind,:),'ok')
hold on
scatter(OSI_all(pref_dir_ind,:), plaidSI_all(pref_dir_ind,:),'or')
ylim([-1 1])
xlim([0 1])
ylabel('Suppression index')
xlabel('OSI')
print(fullfile(summaryDir, 'randPhase_OSIorPrefIndexScatters_8phase.pdf'),'-dpdf','-fillpage')



figure;
scatter(OSI_all(ind,:),testPI_all(ind,:),'ok')
hold on
scatter(OSI_all(pref_dir_ind,:),testPI_all(pref_dir_ind,:),'or')
xlabel('OSI')
ylabel('Stim preference index')
refline(1)
