close all; clear all; clc;
doRedChannel = 0;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'PhaseRevSummary');
ds = 'CrossOriRandPhase_lowSF_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);   
resp_ind_only_all = [];
plaidSI_all = [];
testPI_all = [];
f2overf1_all = [];
totCells = 0;
for iexp = 1:nexp
    if ~isempty(expt(iexp).prFolder)
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        area = expt(iexp).img_loc{1};
        ImgFolder = expt(iexp).coFolder;
        nrun = length(ImgFolder);
        co_run_str = catRunName(cell2mat(ImgFolder), nrun);

        ImgFolder = expt(iexp).prFolder;
        nrun = length(ImgFolder);
        pr_run_str = catRunName(cell2mat(ImgFolder), nrun);
        
        fprintf([mouse ' ' date '\n'])
        
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_respData.mat']))
        load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' pr_run_str], [date '_' mouse '_' pr_run_str '_f1f2.mat']))
        
        
        resp_ind_only_all = [resp_ind_only_all; intersect(resp_ind,resp_ind_phase)+totCells];
        plaid_resp = mean(resp_cell{end,end,1},2);
        mask_resp = mean(resp_cell{end,1,1},2);
        test_resp = mean(resp_cell{1,end,1},2);
        plaid_resp(find(plaid_resp<0)) = 0;
        mask_resp(find(mask_resp<0)) = 0;
        test_resp(find(test_resp<0)) = 0;
        
        plaidSI_all = [plaidSI_all; (plaid_resp-(mask_resp+test_resp)) ./ (plaid_resp + mask_resp + test_resp)];
        testPI_all = [testPI_all; abs((test_resp-mask_resp) ./ (mask_resp+test_resp))];
        f2overf1_all = [f2overf1_all f2overf1];    
        totCells = totCells+size(plaid_resp,1);
    end
end

figure; 
subplot(2,2,1)
scatter(plaidSI_all(resp_ind_only_all),f2overf1_all(resp_ind_only_all))
xlim([-1 1])
ylim([0 2.5])
xlabel('Suppression index')
ylabel('F2/F1')
subplot(2,2,2)
scatter(plaidSI_all(resp_ind_only_all),testPI_all(resp_ind_only_all))
xlim([-1 1])
ylim([0 1])
xlabel('Suppression index')
ylabel('Stim Pref index')
subplot(2,2,3)
scatter(testPI_all(resp_ind_only_all),f2overf1_all(resp_ind_only_all))
xlim([0 1])
ylim([0 2.5])
xlabel('Stim Pref index')
ylabel('F2/F1')
suptitle([date ' ' mouse '- Phase reversal'])
print(fullfile(summaryDir, 'phaseRevSummary_f2f1_SI_PI_scatters.pdf'),'-dpdf','-bestfit');


figure;
subplot(2,2,1)
cdfplot(f2overf1_all(intersect(resp_ind_only_all, find(plaidSI_all<0))))
hold on
cdfplot(f2overf1_all(intersect(resp_ind_only_all, find(plaidSI_all>0))))
xlabel('F2/F1')
title('')
legend({'SI<0','SI>0'})
subplot(2,2,2)
cdfplot(f2overf1_all(intersect(resp_ind_only_all, find(testPI_all<0.5))))
hold on
cdfplot(f2overf1_all(intersect(resp_ind_only_all, find(testPI_all>0.5))))
hold on
xlabel('F2/F1')
title('')
legend({'PI<0.5','PI>0.5'})
subplot(2,2,3)
cdfplot(testPI_all(intersect(resp_ind_only_all, find(plaidSI_all<0))))
hold on
cdfplot(testPI_all(intersect(resp_ind_only_all, find(plaidSI_all>0))))
hold on
xlabel('StimPref')
title('')
legend({'SI<0','SI>0'})
suptitle(['Phase reversal summary- n = ' num2str(length(resp_ind_only_all))])
print(fullfile(summaryDir, 'phaseRevSummary_f2f1_SI_PI_cdfs.pdf'),'-dpdf','-bestfit');

figure;
subplot(2,2,1)
d = 1;
cdfplot(plaidSI_all(intersect(resp_ind_only_all, find(f2overf1_all<d))))
hold on
cdfplot(plaidSI_all(intersect(resp_ind_only_all, find(f2overf1_all>d))))
xlabel('SI')
title('')
legend({['F2/F1<' num2str(d)],['F2/F1>' num2str(d)]},'location','southeast')
subplot(2,2,2)
cdfplot(testPI_all(intersect(resp_ind_only_all, find(f2overf1_all<d))))
hold on
cdfplot(testPI_all(intersect(resp_ind_only_all, find(f2overf1_all>d))))
hold on
xlabel('PI')
title('')
legend({['F2/F1<' num2str(d)],['F2/F1>' num2str(d)]},'location','southeast')
subplot(2,2,3)
d = 0.4;
cdfplot(plaidSI_all(intersect(resp_ind_only_all, find(f2overf1_all<d))))
hold on
cdfplot(plaidSI_all(intersect(resp_ind_only_all, find(f2overf1_all>d))))
xlabel('SI')
title('')
legend({['F2/F1<' num2str(d)],['F2/F1>' num2str(d)]},'location','southeast')
subplot(2,2,4)
cdfplot(testPI_all(intersect(resp_ind_only_all, find(f2overf1_all<d))))
hold on
cdfplot(testPI_all(intersect(resp_ind_only_all, find(f2overf1_all>d))))
hold on
xlabel('PI')
title('')
legend({['F2/F1<' num2str(d)],['F2/F1>' num2str(d)]},'location','southeast')

suptitle(['Phase reversal summary- n = ' num2str(length(resp_ind_only_all))])
print(fullfile(summaryDir, 'phaseRevSummary_f2f1_cdfs.pdf'),'-dpdf','-bestfit');

edges = [-1.1:.2:1.1 ];
plaidSI_all_resp = plaidSI_all(resp_ind_only_all);
f1f2_all_resp = f2overf1_all(resp_ind_only_all)';
[n edges bin] = histcounts(plaidSI_all_resp,edges);
f1f2_avg = zeros(length(n),2);
plaidSI_avg = zeros(length(n),2);
for ibin = 1:length(n)
    ind = find(bin==ibin);
    f1f2_avg(ibin,1) = median(f1f2_all_resp(ind,1),1);
    f1f2_avg(ibin,2) = std(f1f2_all_resp(ind,1),[],1);
    plaidSI_avg(ibin,1) = median(plaidSI_all_resp(ind,1),1);
    plaidSI_avg(ibin,2) = std(plaidSI_all_resp(ind,1), [],1);
end
figure; 
errorbar(plaidSI_avg(:,1),f1f2_avg(:,1), f1f2_avg(:,2),f1f2_avg(:,2),plaidSI_avg(:,2),plaidSI_avg(:,2),'o')
ylim([0 1])
