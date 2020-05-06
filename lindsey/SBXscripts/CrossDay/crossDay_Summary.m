crossDay_16Dir_SLC
prefOri_D1_all = [];
prefOri_D2_all = [];
maxResp_D1_all = [];
maxResp_D2_all = [];
prefOri_diff_all = [];
maxResp_diff_all = [];
goodfit_all = [];
D1_fitRel_all = [];
mice = [];
for iexp = 1:size(expt,2)
    ireg = 1;
    mouse = expt(iexp).mouse;
    reg_date = expt(iexp).reg_dates(ireg,:);
    reg_run = expt(iexp).reg_runs(ireg,:);
    frame_rate = 15.5;
    run_str = catRunName(reg_run, 1);
    fnout = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P';
    load(fullfile(fnout, [reg_date '_' mouse], [reg_date '_' mouse '_' run_str], [reg_date '_' mouse '_' run_str '_compAcrossDaysTuning.mat']));
    goodfit_all = [goodfit_all goodfit+size(prefOri_D1_all,2)];
    prefOri_D1_all = [prefOri_D1_all prefOri_D1];
    prefOri_D2_all = [prefOri_D2_all prefOri_D2];
    maxResp_D1_all = [maxResp_D1_all maxResp_D1];
    maxResp_D2_all = [maxResp_D2_all maxResp_D2];
    prefOri_diff_all = [prefOri_diff_all prefOri_diff];
    maxResp_diff_all = [maxResp_diff_all maxResp_diff];
    D1_fitRel_all = [D1_fitRel_all D1_fitReliability];
    mice{iexp} = mouse;
end

figure; 
subplot(2,3,1)
scatter(prefOri_D1_all(goodfit_all),prefOri_D2_all(goodfit_all));
axis square
xlim([0 180])
ylim([0 180])
refline(1,0)
xlabel('Day 1 Pref Ori')
ylabel('Day 2 Pref Ori')

subplot(2,3,2)
bins = [5:10:95];
hist(prefOri_diff_all,bins)
xlim([0 90])
axis square
xlabel('D1 vs D2 Pref ori')
ylabel('Number of cells')

subplot(2,3,3)
hist(D1_fitRel_all(goodfit_all),bins)
xlim([0 90])
axis square
xlabel('D1 90% CI')
ylabel('Number of cells')

subplot(2,3,4)
scatter(maxResp_D1_all(goodfit_all),prefOri_diff_all)
xlabel('Day 1 max dF/F')
ylabel('D1 vs D2 Pref ori')
axis square

subplot(2,3,5)
scatter(maxResp_diff_all,prefOri_diff_all)
xlabel('D1 vs D2 max dF/F')
ylabel('D1 vs D2 Pref ori')
axis square

suptitle([cell2mat(mice) '- n = ' num2str(length(goodfit_all)) ' cells'])
fn_summary = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\CrossDay';
print(fullfile(fn_summary, [cell2mat(mice) '_regD' num2str(ireg) '_summaryAcrossDaysTuning.pdf']),'-dpdf', '-bestfit')
