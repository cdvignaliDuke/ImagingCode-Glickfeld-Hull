close all;clear all;clc;

ds = 'CrossOriSingleStimAdapt_ExptList';
rc = behavConstsAV;
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

eval(ds)
nexp = length(expt);

totCells = 0;
preftestonly_ind_all = [];
prefmaskonly_ind_all = [];
prefplaidonly_ind_all = [];
resptest_ind_all = [];
respmask_ind_all = [];
respplaid_ind_all = [];
noadapt_resp_cell_all= cell(5,5);
singadapt_resp_cell_all= cell(5,5);
dir_resp_all = [];
dirtuned_ind_all = [];
dirresp_ind_all = [];
u1hat_all = [];
k1hat_all = [];
R1hat_all = [];
R2hat_all = [];
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    run_str = ['runs-' cell2mat(expt(iexp).coFolder)];
    dir_run_str = ['runs-' cell2mat(expt(iexp).dirFolder)];
    
    fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);
    dir_fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str]);
    load(fullfile(fn,[date '_' mouse '_' run_str '_allCellResp.mat']));
    %load(fullfile(fn,[date '_' mouse '_' run_str '_allCellRespLate.mat']));
    load(fullfile(fn,[date '_' mouse '_' run_str '_dataStim.mat']));
    nMask = length(maskCons);
    nTest = length(testCons);

    for im = 1:nMask
        for it = 1:nTest
            noadapt_resp_cell_all{im,it} = [noadapt_resp_cell_all{im,it}; nanmean(noadapt_resp_cell{im,it},2)];
            singadapt_resp_cell_all{im,it} = [singadapt_resp_cell_all{im,it}; nanmean(singadapt_resp_cell{im,it},2)];
        end
    end
    preftestonly_ind_all = [preftestonly_ind_all; preftestonly_ind+totCells];
    prefmaskonly_ind_all = [prefmaskonly_ind_all; prefmaskonly_ind+totCells];
    prefplaidonly_ind_all = [prefplaidonly_ind_all; prefplaidonly_ind+totCells];
    resptest_ind_all = [resptest_ind_all; resptest_ind+totCells];
    respmask_ind_all = [respmask_ind_all; respmask_ind+totCells];
    respplaid_ind_all = [respplaid_ind_all; respplaid_ind+totCells];
    
    load(fullfile(dir_fn,[date '_' mouse '_' dir_run_str '_stimData.mat']));
    load(fullfile(dir_fn,[date '_' mouse '_' dir_run_str '_dfofData.mat']));
    dir_resp_all = cat(2, dir_resp_all, data_dfof_dir);
    dirresp_ind_all = [dirresp_ind_all  dirresp_ind+totCells];
    dirtuned_ind_all = [dirtuned_ind_all  dirtuned_ind+totCells];
    u1hat_all = [u1hat_all rad2deg(u1_hat)];
    k1hat_all = [k1hat_all k1_hat];
    R1hat_all = [R1hat_all R1_hat];
    R2hat_all = [R2hat_all R2_hat];
    
    totCells = totCells+size(noadapt_resp_cell{1,1},1); 
        
    fprintf([num2str(length(preftestonly_ind)) ' ' num2str(length(prefmaskonly_ind)) ' ' num2str(length(prefplaidonly_ind)) '\n'])
end

prefTestInd = intersect(preftestonly_ind_all,resptest_ind_all);
prefMaskInd = intersect(prefmaskonly_ind_all,respmask_ind_all);
prefPlaidInd = intersect(prefplaidonly_ind_all,respplaid_ind_all);

%% cross ori measures
fnout = fullfile(LG_base, 'Analysis\2P\CrossOri');
if ~exist(fnout)
    mkdir(fnout);
end

mask_only_resp = noadapt_resp_cell_all{5,1}-noadapt_resp_cell_all{1,1};
mask_only_resp(find(mask_only_resp<0)) = 0;
test_only_resp = noadapt_resp_cell_all{1,5}-noadapt_resp_cell_all{1,1};
test_only_resp(find(test_only_resp<0)) = 0;
plaid_only_resp =noadapt_resp_cell_all{5,5}-noadapt_resp_cell_all{1,1};
plaid_only_resp(find(plaid_only_resp<0)) = 0;
MTpref_all = abs((mask_only_resp-test_only_resp)./(mask_only_resp+test_only_resp));
SI_all = (plaid_only_resp-(mask_only_resp+test_only_resp))./(plaid_only_resp+(mask_only_resp+test_only_resp));
resp_all = unique([resptest_ind_all; respmask_ind_all; respplaid_ind_all]);
figure;
edges_SI = [-1:0.2:1];
subplot(2,4,1)
hist(SI_all(resp_all),edges_SI);
vline(mean(SI_all(resp_all)))
title('All responsive')
xlabel('Suppression index')
ylabel('Number of cells')
subplot(2,4,2)
hist(SI_all(resptest_ind_all),edges_SI);
vline(mean(SI_all(resptest_ind_all)))
title('Test responsive')
xlabel('Suppression index')
ylabel('Number of cells')
subplot(2,4,3)
hist(SI_all(respmask_ind_all),edges_SI);
vline(mean(SI_all(respmask_ind_all)))
title('Mask responsive')
xlabel('Suppression index')
ylabel('Number of cells')
subplot(2,4,4)
hist(SI_all(respplaid_ind_all),edges_SI);
vline(mean(SI_all(respplaid_ind_all)))
title('Plaid responsive')
xlabel('Suppression index')
ylabel('Number of cells')
subplot(2,4,6)
hist(SI_all(prefTestInd),edges_SI);
vline(mean(SI_all(prefTestInd)))
title('Pref Test')
xlabel('Suppression index')
ylabel('Number of cells')
subplot(2,4,7)
hist(SI_all(prefMaskInd),edges_SI);
vline(mean(SI_all(prefMaskInd)))
title('Pref Mask')
xlabel('Suppression index')
ylabel('Number of cells')
subplot(2,4,8)
hist(SI_all(prefPlaidInd),edges_SI);
vline(mean(SI_all(prefPlaidInd)))
title('Pref Plaid')
xlabel('Suppression index')
ylabel('Number of cells')
print(fullfile(fnout, 'crossOriSummary_suppIndex.pdf'),'-dpdf','-bestfit');

figure;
subplot(1,2,1)
scatter(MTpref_all(resp_all),SI_all(resp_all))
ylabel('Suppression index')
xlabel('abs((M-T)/(M+T))')
edges_OSI = [0:0.1:1];
[n bin] = histc(MTpref_all,edges_OSI);
subplot(1,2,2)
for i =1:length(n)
    ind = intersect(resp_all,find(bin == i));
    errorbar(edges_OSI(i), mean(SI_all(ind,:),1), std(SI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
end
ylabel('Suppression index')
xlabel('abs((M-T)/(M+T))')
ylim([-1 1])
xlim([-0.05 1.05])
print(fullfile(fnout, 'crossOriSummary_suppIndexVsMTpref.pdf'),'-dpdf','-bestfit');

dir_resp_all_avg = squeeze(mean(dir_resp_all(resp_win,:,:),1));
figure;
subplot(2,3,1)
errorbar(dirs, mean(dir_resp_all_avg(prefTestInd,:),1), std(dir_resp_all_avg(prefTestInd,:),[],1)./sqrt(length(prefTestInd)),'ok')
fprintf(['Pref mask max = ' num2str(dirs(find(mean(dir_resp_all_avg(prefTestInd,:),1) == max(mean(dir_resp_all_avg(prefTestInd,:),1),[],2)))) '\n'])
title('Pref Test')
xlabel('Direction')
ylabel('dF/F')
subplot(2,3,2)
errorbar(dirs, mean(dir_resp_all_avg(prefMaskInd,:),1), std(dir_resp_all_avg(prefMaskInd,:),[],1)./sqrt(length(prefMaskInd)),'ok')
fprintf(['Pref mask max = ' num2str(dirs(find(mean(dir_resp_all_avg(prefMaskInd,:),1) == max(mean(dir_resp_all_avg(prefMaskInd,:),1),[],2)))) '\n'])
title('Pref Mask')
xlabel('Direction')
ylabel('dF/F')
subplot(2,3,3)
errorbar(dirs, mean(dir_resp_all_avg(prefPlaidInd,:),1), std(dir_resp_all_avg(prefPlaidInd,:),[],1)./sqrt(length(prefPlaidInd)),'ok')
fprintf(['Pref mask max = ' num2str(dirs(find(mean(dir_resp_all_avg(prefPlaidInd,:),1) == max(mean(dir_resp_all_avg(prefPlaidInd,:),1),[],2)))) '\n'])
title('Pref Plaid')
xlabel('Direction')
ylabel('dF/F')
subplot(2,3,4)
hist(u1hat_all(prefTestInd))
xlabel('Preferred direction')
ylabel('Number of cells')
subplot(2,3,5)
hist(u1hat_all(prefMaskInd))
xlabel('Preferred direction')
ylabel('Number of cells')
subplot(2,3,6)
hist(u1hat_all(prefPlaidInd))
xlabel('Preferred direction')
ylabel('Number of cells')
print(fullfile(fnout, 'crossOriSummary_distPrefDirBySI.pdf'),'-dpdf','-bestfit');

[max_val max_ind] = max(dir_resp_all_avg,[],2);
max_dir = dirs(max_ind);
DSI_fit = (R1hat_all-R2hat_all)./(R1hat_all+R2hat_all);
nullDir = max_ind+nDir/2;
nullDir(find(nullDir>nDir)) = nullDir(find(nullDir>nDir))-nDir;
min_val = zeros(totCells,1);
for iCell = 1:totCells
    min_val(iCell,1) = dir_resp_all_avg(iCell,nullDir(iCell));
end
min_val(find(min_val<0)) = 0;
max_val(find(max_val<0)) = 0;
DSI_resp = (max_val-min_val)./(max_val+min_val);
ori_resp_all_avg = mean(reshape(dir_resp_all_avg,[totCells nDir/2 2]),3);
oris = dirs(1:nDir/2);
[max_val max_ind] = max(ori_resp_all_avg,[],2);
nullOri = max_ind+nDir/4;
nullOri(find(nullOri>nDir/2)) = nullOri(find(nullOri>nDir/2))-nDir/2;
min_val = zeros(totCells,1);
for iCell = 1:totCells
    min_val(iCell,1) = ori_resp_all_avg(iCell,nullOri(iCell));
end
min_val(find(min_val<0)) = 0;
max_val(find(max_val<0)) = 0;
OSI_resp = (max_val-min_val)./(max_val+min_val);

figure;
subplot(3,2,1)
scatter(DSI_fit(resp_all),SI_all(resp_all))
ylabel('Suppression index')
xlabel('DSI-fit')
edges_DSI = [0:0.1:1];
[n bin] = histc(DSI_fit,edges_DSI);
subplot(3,2,2)
for i =1:length(n)
    ind = intersect(resp_all,find(bin == i));
    errorbar(edges_DSI(i), mean(SI_all(ind,:),1), std(SI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
end
ylabel('Suppression index')
xlabel('DSI-fit')
ylim([-1 1])
xlim([-0.05 1.05])
subplot(3,2,3)
scatter(DSI_resp(resp_all),SI_all(resp_all))
ylabel('Suppression index')
xlabel('DSI-resp')
edges_DSI = [0:0.1:1];
[n bin] = histc(DSI_resp,edges_DSI);
subplot(3,2,4)
for i =1:length(n)
    ind = intersect(resp_all,find(bin == i));
    errorbar(edges_DSI(i), mean(SI_all(ind,:),1), std(SI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
end
ylabel('Suppression index')
xlabel('DSI-resp')
ylim([-1 1])
xlim([-0.05 1.05])
subplot(3,2,5)
scatter(OSI_resp(resp_all),SI_all(resp_all))
ylabel('Suppression index')
xlabel('OSI-resp')
edges_OSI = [0:0.1:1];
[n bin] = histc(OSI_resp,edges_DSI);
subplot(3,2,6)
for i =1:length(n)
    ind = intersect(resp_all,find(bin == i));
    errorbar(edges_OSI(i), mean(SI_all(ind,:),1), std(SI_all(ind,:),[],1)./sqrt(length(ind)),'ok')
    hold on
end
ylabel('Suppression index')
xlabel('OSI-resp')
ylim([-1 1])
xlim([-0.05 1.05])
print(fullfile(fnout, 'crossOriSummary_suppIndexVsDSIorOSI.pdf'),'-dpdf','-bestfit');

figure;
subplot(2,3,1)
cdfplot(DSI_resp(resptest_ind_all))
hold on
cdfplot(DSI_resp(respmask_ind_all))
cdfplot(DSI_resp(respplaid_ind_all))
xlabel('DSI')
title('')
legend({'Resp Test', 'Resp Mask', 'Resp Plaid'}, 'location', 'southeast');
subplot(2,3,2)
cdfplot(OSI_resp(resptest_ind_all))
hold on
cdfplot(OSI_resp(respmask_ind_all))
cdfplot(OSI_resp(respplaid_ind_all))
xlabel('OSI')
title('')
subplot(2,3,3)
cdfplot(MTpref_all(resptest_ind_all))
hold on
cdfplot(MTpref_all(respmask_ind_all))
cdfplot(MTpref_all(respplaid_ind_all))
xlabel('M-T/M+T')
title('')
subplot(2,3,4)
cdfplot(DSI_resp(prefTestInd))
hold on
cdfplot(DSI_resp(prefMaskInd))
cdfplot(DSI_resp(prefPlaidInd))
xlabel('DSI')
title('')
legend({'Pref Test', 'Pref Mask', 'Pref Plaid'}, 'location', 'southeast');
subplot(2,3,5)
cdfplot(OSI_resp(prefTestInd))
hold on
cdfplot(OSI_resp(prefMaskInd))
cdfplot(OSI_resp(prefPlaidInd))
xlabel('OSI')
title('')
subplot(2,3,6)
cdfplot(MTpref_all(prefTestInd))
hold on
cdfplot(MTpref_all(prefMaskInd))
cdfplot(MTpref_all(prefPlaidInd))
xlabel('M-T/M+T')
title('')
print(fullfile(fnout, 'crossOriSummary_selectivityByPref&Resp.pdf'),'-dpdf','-bestfit');

figure;
subplot(2,2,1)
scatter(SI_all(resptest_ind_all),u1hat_all(resptest_ind_all))
ylabel('Pref Dir from fit')
xlabel('Suppression index')
title('Test')
subplot(2,2,2)
scatter(SI_all(respmask_ind_all),u1hat_all(respmask_ind_all))
ylabel('Pref Dir from fit')
xlabel('Suppression index')
title('Mask')
subplot(2,2,3)
scatter(SI_all(respplaid_ind_all),u1hat_all(respplaid_ind_all))
ylabel('Pref Dir from fit')
xlabel('Suppression index')
title('Plaid')
subplot(2,2,4)
scatter(SI_all(prefTestInd),u1hat_all(prefTestInd))
hold on
scatter(SI_all(prefMaskInd),u1hat_all(prefMaskInd))
scatter(SI_all(prefPlaidInd),u1hat_all(prefPlaidInd))
ylabel('Pref Dir from fit')
xlabel('Suppression index')
print(fullfile(fnout, 'crossOriSummary_prefDirBySI_allPref&Resp.pdf'),'-dpdf','-bestfit');

%%
prefTestOnly_noadapt_resp = zeros(nMask,nTest,2);
prefMaskOnly_noadapt_resp = zeros(nMask,nTest,2);
prefPlaidOnly_noadapt_resp = zeros(nMask,nTest,2);
prefTestOnly_singadapt_resp = zeros(nMask,nTest,2);
prefMaskOnly_singadapt_resp = zeros(nMask,nTest,2);
prefPlaidOnly_singadapt_resp = zeros(nMask,nTest,2);
respMask_noadapt_resp = zeros(nMask,nTest,2);
respMask_singadapt_resp = zeros(nMask,nTest,2);

for im = 1:nMask
    for it = 1:nTest
        prefTestOnly_noadapt_resp(im,it,1) = nanmean(noadapt_resp_cell_all{im,it}(prefTestInd,:),1);
        prefTestOnly_noadapt_resp(im,it,2) = nanstd(noadapt_resp_cell_all{im,it}(prefTestInd,:),[],1)./sqrt(length(prefTestInd));
        prefMaskOnly_noadapt_resp(im,it,1) = nanmean(noadapt_resp_cell_all{im,it}(prefMaskInd,:),1);
        prefMaskOnly_noadapt_resp(im,it,2) = nanstd(noadapt_resp_cell_all{im,it}(prefMaskInd,:),[],1)./sqrt(length(prefMaskInd));
        prefPlaidOnly_noadapt_resp(im,it,1) = nanmean(noadapt_resp_cell_all{im,it}(prefPlaidInd,:),1);
        prefPlaidOnly_noadapt_resp(im,it,2) = nanstd(noadapt_resp_cell_all{im,it}(prefPlaidInd,:),[],1)./sqrt(length(prefPlaidInd));
        prefTestOnly_singadapt_resp(im,it,1) = nanmean(singadapt_resp_cell_all{im,it}(prefTestInd,:),1);
        prefTestOnly_singadapt_resp(im,it,2) = nanstd(singadapt_resp_cell_all{im,it}(prefTestInd,:),[],1)./sqrt(length(prefTestInd));
        prefMaskOnly_singadapt_resp(im,it,1) = nanmean(singadapt_resp_cell_all{im,it}(prefMaskInd,:),1);
        prefMaskOnly_singadapt_resp(im,it,2) = nanstd(singadapt_resp_cell_all{im,it}(prefMaskInd,:),[],1)./sqrt(length(prefMaskInd));
        prefPlaidOnly_singadapt_resp(im,it,1) = nanmean(singadapt_resp_cell_all{im,it}(prefPlaidInd,:),1);
        prefPlaidOnly_singadapt_resp(im,it,2) = nanstd(singadapt_resp_cell_all{im,it}(prefPlaidInd,:),[],1)./sqrt(length(prefPlaidInd));
        respMask_noadapt_resp(im,it,1) = nanmean(noadapt_resp_cell_all{im,it}(respmask_ind_all,:),1);
        respMask_noadapt_resp(im,it,2) = nanstd(noadapt_resp_cell_all{im,it}(respmask_ind_all,:),[],1)./sqrt(length(respmask_ind_all));
        respMask_singadapt_resp(im,it,1) = nanmean(singadapt_resp_cell_all{im,it}(respmask_ind_all,:),1);
        respMask_singadapt_resp(im,it,2) = nanstd(singadapt_resp_cell_all{im,it}(respmask_ind_all,:),[],1)./sqrt(length(respmask_ind_all));
    end
end

% TestOnly
figure;
subplot(1,2,1)
for im = 1:nMask
    errorbar(testCons, prefTestOnly_noadapt_resp(im,:,1),prefTestOnly_noadapt_resp(im,:,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,2,2)
for im = 1:nMask
    errorbar(testCons, prefTestOnly_singadapt_resp(im,:,1),prefTestOnly_singadapt_resp(im,:,2),'-o')
    hold on
end
title('SingleStim')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend(num2str(maskCons'))
suptitle(['Test preferring cells- n = ' num2str(length(prefTestInd))])
print(fullfile(fnout, 'adaptSummary_testPrefAll.pdf'),'-dpdf','-bestfit');

%Mask Only
figure;
subplot(1,2,1)
for it = 1:nTest
    errorbar(maskCons, prefMaskOnly_noadapt_resp(:,it,1),prefMaskOnly_noadapt_resp(:,it,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,2,2)
for it = 1:nTest
    errorbar(maskCons, prefMaskOnly_singadapt_resp(:,it,1),prefMaskOnly_singadapt_resp(:,it,2),'-o')
    hold on
end
title('SingleStim')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend(num2str(testCons'))
suptitle(['Mask preferring cells- n = ' num2str(length(prefMaskInd))])
print(fullfile(fnout, 'adaptSummary_maskPrefAll.pdf'),'-dpdf','-bestfit');

% Plaid Only
figure;
subplot(1,2,1)
for im = 1:nMask
    errorbar(testCons, prefPlaidOnly_noadapt_resp(im,:,1),prefPlaidOnly_noadapt_resp(im,:,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,2,2)
for im = 1:nMask
    errorbar(testCons, prefPlaidOnly_singadapt_resp(im,:,1),prefPlaidOnly_singadapt_resp(im,:,2),'-o')
    hold on
end
title('SingleStim')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend(num2str(maskCons'))
suptitle(['Plaid preferring cells- n = ' num2str(length(prefPlaidInd))])
print(fullfile(fnout, 'adaptSummary_plaidPrefOnly.pdf'),'-dpdf','-bestfit');

%Mask responsive
figure;
subplot(1,2,1)
for it = 1:nTest
    errorbar(testCons, respMask_noadapt_resp(:,it,1),respMask_noadapt_resp(:,it,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,2,2)
for it = 1:nTest
    errorbar(testCons, respMask_singadapt_resp(:,it,1),respMask_singadapt_resp(:,it,2),'-o')
    hold on
end
title('SingleStim')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend(num2str(maskCons'))
suptitle(['Mask responsive cells- n = ' num2str(length(respmask_ind_all))])
print(fullfile(fnout, 'adaptSummary_respMaskAll.pdf'),'-dpdf','-bestfit');

%by mask
figure;
for it = 1:nTest
    subplot(2,3,it)
    errorbar(maskCons, respMask_noadapt_resp(:,it,1),respMask_noadapt_resp(:,it,2),'-o')
    hold on
    errorbar(maskCons, respMask_singadapt_resp(:,it,1),respMask_singadapt_resp(:,it,2),'-o')
    title(['Test = ' num2str(testCons(it))])
    ylabel('dF/F')
    xlabel('Mask Contrast')
    ylim([-0.1 0.5])
end
legend('No adapt', 'SingleStim');
suptitle(['Mask responsive cells- n = ' num2str(length(respmask_ind_all))])
print(fullfile(fnout, 'adaptSummary_byMask_respMaskAll.pdf'),'-dpdf','-bestfit');

figure;
for im = 1:nMask
    subplot(2,3,im)
    errorbar(testCons, prefTestOnly_noadapt_resp(im,:,1),prefTestOnly_noadapt_resp(im,:,2),'-o')
    hold on
    errorbar(testCons, prefTestOnly_singadapt_resp(im,:,1),prefTestOnly_singadapt_resp(im,:,2),'-o')
    title(['Mask = ' num2str(maskCons(im))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
end
legend('No adapt', 'SingleStim');
suptitle(['Test only preferring cells- n = ' num2str(length([prefTestInd]))])
print(fullfile(fnout, 'adaptSummary_byMask_prefTestAll.pdf'),'-dpdf','-bestfit');


figure;
for it = 1:nTest
    subplot(2,3,it)
    errorbar(maskCons, prefMaskOnly_noadapt_resp(:,it,1),prefMaskOnly_noadapt_resp(:,it,2),'-o')
    hold on
    errorbar(maskCons, prefMaskOnly_singadapt_resp(:,it,1),prefMaskOnly_singadapt_resp(:,it,2),'-o')
    title(['Test = ' num2str(testCons(it))])
    ylabel('dF/F')
    xlabel('Mask Contrast')
    ylim([-0.1 0.5])
end
legend('No adapt', 'SingleStim');
suptitle(['Mask only preferring cells- n = ' num2str(length([prefMaskInd]))])
print(fullfile(fnout, 'adaptSummary_byMask_prefMaskAll.pdf'),'-dpdf','-bestfit');

figure;
for im = 1:nMask
    subplot(2,3,im)
    errorbar(testCons, prefPlaidOnly_noadapt_resp(im,:,1),prefPlaidOnly_noadapt_resp(im,:,2),'-o')
    hold on
    errorbar(testCons, prefPlaidOnly_singadapt_resp(im,:,1),prefPlaidOnly_singadapt_resp(im,:,2),'-o')
    title(['Mask = ' num2str(maskCons(im))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
end
legend('No adapt', 'SingleStim');
suptitle(['Plaid only preferring cells- n = ' num2str(length([prefPlaidInd]))])
print(fullfile(fnout, 'adaptSummary_byMask_prefPlaidOnly.pdf'),'-dpdf','-bestfit');

%%
prefTestOnly_noadapt_respCells = zeros(length(prefTestInd),nMask,nTest);
prefMaskOnly_noadapt_respCells = zeros(length(prefMaskInd),nMask,nTest);
prefPlaidOnly_noadapt_respCells = zeros(length(prefPlaidInd),nMask,nTest);
respMask_noadapt_respCells = zeros(length(respmask_ind_all),nMask,nTest);
prefTestOnly_singadapt_respCells = zeros(length(prefTestInd),nMask,nTest);
prefMaskOnly_singadapt_respCells = zeros(length(prefMaskInd),nMask,nTest);
prefPlaidOnly_singadapt_respCells = zeros(length(prefPlaidInd),nMask,nTest);
respMask_singadapt_respCells = zeros(length(respmask_ind_all),nMask,nTest);
for im = 1:nMask
    for it = 1:nTest
        prefTestOnly_noadapt_respCells(:,im,it) = noadapt_resp_cell_all{im,it}(prefTestInd,:)-noadapt_resp_cell_all{1,1}(prefTestInd,:);
        prefMaskOnly_noadapt_respCells(:,im,it) = noadapt_resp_cell_all{it,im}(prefMaskInd,:)-noadapt_resp_cell_all{1,1}(prefMaskInd,:);
        prefPlaidOnly_noadapt_respCells(:,im,it) = noadapt_resp_cell_all{im,it}(prefPlaidInd,:)-noadapt_resp_cell_all{1,1}(prefPlaidInd,:);
        prefTestOnly_singadapt_respCells(:,im,it) = singadapt_resp_cell_all{im,it}(prefTestInd,:)-singadapt_resp_cell_all{1,1}(prefTestInd,:);
        prefMaskOnly_singadapt_respCells(:,im,it) = singadapt_resp_cell_all{it,im}(prefMaskInd,:)-singadapt_resp_cell_all{1,1}(prefMaskInd,:);
        prefPlaidOnly_singadapt_respCells(:,im,it) = singadapt_resp_cell_all{im,it}(prefPlaidInd,:)-singadapt_resp_cell_all{1,1}(prefPlaidInd,:);
        respMask_noadapt_respCells(:,im,it) = noadapt_resp_cell_all{it,im}(respmask_ind_all,:)-noadapt_resp_cell_all{1,1}(respmask_ind_all,:);
        respMask_singadapt_respCells(:,im,it) = singadapt_resp_cell_all{it,im}(respmask_ind_all,:)-singadapt_resp_cell_all{1,1}(respmask_ind_all,:); 
    end
end

%set negative values to 0
prefTestOnly_noadapt_respCells(find(prefTestOnly_noadapt_respCells<0)) = 0;
prefMaskOnly_noadapt_respCells(find(prefMaskOnly_noadapt_respCells<0)) = 0;
prefPlaidOnly_noadapt_respCells(find(prefPlaidOnly_noadapt_respCells<0)) = 0;
respMask_noadapt_respCells(find(respMask_noadapt_respCells<0)) = 0;
prefTestOnly_singadapt_respCells(find(prefTestOnly_singadapt_respCells<0)) = 0;
prefMaskOnly_singadapt_respCells(find(prefMaskOnly_singadapt_respCells<0)) = 0;
prefPlaidOnly_singadapt_respCells(find(prefPlaidOnly_singadapt_respCells<0)) = 0;
respMask_singadapt_respCells(find(respMask_singadapt_respCells<0)) = 0;

figure;
subplot(2,2,1)
SI = (prefTestOnly_singadapt_respCells-prefTestOnly_noadapt_respCells)./(prefTestOnly_noadapt_respCells+prefTestOnly_singadapt_respCells);
imagesc(squeeze(nanmean(SI,1)))
colormap redblue
clim([-1 1])
axis square
ylabel('Test Contrast')
set(gca,'Xtick',[1:5])
set(gca,'Ytick',[1:5])
xticklabels(chop(maskCons,2))
xlabel('Mask Contrast')
yticklabels(chop(testCons,2))
colorbar
title({'Test preferring cells- ' , ['n = ' num2str(length(prefTestInd))]})
subplot(2,2,2)
SI = (prefMaskOnly_singadapt_respCells-prefMaskOnly_noadapt_respCells)./(prefMaskOnly_noadapt_respCells+prefMaskOnly_singadapt_respCells);
imagesc(squeeze(nanmean(SI,1)))
colormap redblue
clim([-1 1])
axis square
ylabel('Test Contrast')
set(gca,'Xtick',[1:5])
set(gca,'Ytick',[1:5])
xticklabels(chop(maskCons,2))
xlabel('Mask Contrast')
yticklabels(chop(testCons,2))
colorbar
title({'Mask preferring cells- ' , ['n = ' num2str(length(prefMaskInd))]})
subplot(2,2,3)
SI = (prefPlaidOnly_singadapt_respCells-prefPlaidOnly_noadapt_respCells)./(prefPlaidOnly_noadapt_respCells+prefPlaidOnly_singadapt_respCells);
imagesc(squeeze(nanmean(SI,1)))
colormap redblue
clim([-1 1])
axis square
ylabel('Test Contrast')
set(gca,'Xtick',[1:5])
set(gca,'Ytick',[1:5])
xticklabels(chop(maskCons,2))
xlabel('Mask Contrast')
yticklabels(chop(testCons,2))
colorbar
title({'Plaid preferring cells- ' , ['n = ' num2str(length(prefPlaidInd))]})
subplot(2,2,4)
SI = (respMask_singadapt_respCells-respMask_noadapt_respCells)./(respMask_noadapt_respCells+respMask_singadapt_respCells);
imagesc(squeeze(nanmean(SI,1)))
colormap redblue
clim([-1 1])
axis square
ylabel('Test Contrast')
set(gca,'Xtick',[1:5])
set(gca,'Ytick',[1:5])
xticklabels(chop(maskCons,2))
xlabel('Mask Contrast')
yticklabels(chop(testCons,2))
colorbar
title({'Mask responsive cells- ' , ['n = ' num2str(length(respmask_ind_all))]})
suptitle('Adaptation index')
print(fullfile(fnout, 'adaptSummary_adaptDiff.pdf'),'-dpdf','-bestfit');


%% modulation index
prefTestOnly_noadapt_MI = zeros(length(prefTestInd),nTest,nMask-1);
prefMaskOnly_noadapt_MI = zeros(length(prefMaskInd),nTest,nMask-1);
prefPlaidOnly_noadapt_MI = zeros(length(prefPlaidInd),nTest,nMask-1);
respMask_noadapt_MI = zeros(length(respmask_ind_all),nTest,nMask-1);
prefTestOnly_singadapt_MI = zeros(length(prefTestInd),nTest,nMask-1);
prefMaskOnly_singadapt_MI = zeros(length(prefMaskInd),nTest,nMask-1);
prefPlaidOnly_singadapt_MI = zeros(length(prefPlaidInd),nTest,nMask-1);
respMask_singadapt_MI = zeros(length(respmask_ind_all),nTest,nMask-1);
    
for im = 2:nMask
    prefTestOnly_noadapt_MI(:,:,im-1) = squeeze(-prefTestOnly_noadapt_respCells(:,1,:)+prefTestOnly_noadapt_respCells(:,im,:))./squeeze(prefTestOnly_noadapt_respCells(:,1,:)+prefTestOnly_noadapt_respCells(:,im,:));
    prefMaskOnly_noadapt_MI(:,:,im-1) = squeeze(-prefMaskOnly_noadapt_respCells(:,1,:)+prefMaskOnly_noadapt_respCells(:,im,:))./squeeze(prefMaskOnly_noadapt_respCells(:,1,:)+prefMaskOnly_noadapt_respCells(:,im,:));
    prefPlaidOnly_noadapt_MI(:,:,im-1) = squeeze(-prefPlaidOnly_noadapt_respCells(:,1,:)+prefPlaidOnly_noadapt_respCells(:,im,:))./squeeze(prefPlaidOnly_noadapt_respCells(:,1,:)+prefPlaidOnly_noadapt_respCells(:,im,:));
    respMask_noadapt_MI(:,:,im-1) = squeeze(-respMask_noadapt_respCells(:,1,:)+respMask_noadapt_respCells(:,im,:))./squeeze(respMask_noadapt_respCells(:,1,:)+respMask_noadapt_respCells(:,im,:));    
    
    prefTestOnly_singadapt_MI(:,:,im-1) = squeeze(-prefTestOnly_singadapt_respCells(:,1,:)+prefTestOnly_singadapt_respCells(:,im,:))./squeeze(prefTestOnly_singadapt_respCells(:,1,:)+prefTestOnly_singadapt_respCells(:,im,:));
    prefMaskOnly_singadapt_MI(:,:,im-1) = squeeze(-prefMaskOnly_singadapt_respCells(:,1,:)+prefMaskOnly_singadapt_respCells(:,im,:))./squeeze(prefMaskOnly_singadapt_respCells(:,1,:)+prefMaskOnly_singadapt_respCells(:,im,:));
    prefPlaidOnly_singadapt_MI(:,:,im-1) = squeeze(-prefPlaidOnly_singadapt_respCells(:,1,:)+prefPlaidOnly_singadapt_respCells(:,im,:))./squeeze(prefPlaidOnly_singadapt_respCells(:,1,:)+prefPlaidOnly_singadapt_respCells(:,im,:));
    respMask_singadapt_MI(:,:,im-1) = squeeze(-respMask_singadapt_respCells(:,1,:)+respMask_singadapt_respCells(:,im,:))./squeeze(respMask_singadapt_respCells(:,1,:)+respMask_singadapt_respCells(:,im,:));    
end

figure; 
subplot(2,2,1)
imagesc(squeeze(nanmean(prefTestOnly_noadapt_MI(:,2:5,:),1))); clim([-1 1]); colormap('redblue'); title('Test Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
axis square
subplot(2,2,2)
imagesc(squeeze(nanmean(prefMaskOnly_noadapt_MI(:,2:5,:),1))); clim([-1 1]); colormap('redblue'); title('Mask Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
axis square
subplot(2,2,3)
imagesc(squeeze(nanmean(prefPlaidOnly_noadapt_MI(:,2:5,:),1))); clim([-1 1]); colormap('redblue'); title('Plaid Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
axis square
subplot(2,2,4)
imagesc(squeeze(nanmean(respMask_noadapt_MI(:,2:5,:),1))); clim([-1 1]); colormap('redblue'); title('Mask Resp')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
axis square
print(fullfile(fnout, 'adaptSummary_modulationByPref.pdf'),'-dpdf','-bestfit');

%Test pref
figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefTestOnly_noadapt_MI(:,it,:)),1);
    n_singadapt = sum(~isnan(prefTestOnly_singadapt_MI(:,it,:)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOnly_noadapt_MI(:,it,:),1)), squeeze(nanstd(prefTestOnly_noadapt_MI(:,it,:),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOnly_singadapt_MI(:,it,:),1)), squeeze(nanstd(prefTestOnly_singadapt_MI(:,it,:),1)./sqrt(n_singadapt)), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('MI')
    ylim([-1 1])
    hline(0)
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefTestOnly_noadapt_MI,2)),1);
n_singadapt = sum(~isnan(nanmean(prefTestOnly_singadapt_MI,2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOnly_noadapt_MI,2),1)), squeeze(nanstd(nanmean(prefTestOnly_noadapt_MI,2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOnly_singadapt_MI,2),1)), squeeze(nanstd(nanmean(prefTestOnly_singadapt_MI,2),1)./sqrt(n_singadapt)), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('MI')
ylim([-1 1])
hline(0)
subplot(3,2,6)
n_singadapt = squeeze(sum(~isnan(nanmean(prefTestOnly_singadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefTestOnly_singadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefTestOnly_singadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2),[],1))./sqrt(n_singadapt)),'-o');
title('All Test')
xlabel('Mask contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
hline(0)
legend({'No adapt', 'Single'})
suptitle(['Test Pref- n = ' num2str(min(n_singadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_MI_prefTestAll.pdf'),'-dpdf','-bestfit');
%Mask pref
figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefMaskOnly_noadapt_MI(:,it,:)),1);
    n_singadapt = sum(~isnan(prefMaskOnly_singadapt_MI(:,it,:)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefMaskOnly_noadapt_MI(:,it,:),1)), squeeze(nanstd(prefMaskOnly_noadapt_MI(:,it,:),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefMaskOnly_singadapt_MI(:,it,:),1)), squeeze(nanstd(prefMaskOnly_singadapt_MI(:,it,:),1)./sqrt(n_singadapt)), '-o')
    title(['Mask = ' num2str(testCons(it))])
    xlabel('Test contrast')
    ylabel('MI')
    ylim([-1 1])
    hline(0)
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefMaskOnly_noadapt_MI,2)),1);
n_singadapt = sum(~isnan(nanmean(prefMaskOnly_singadapt_MI,2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefMaskOnly_noadapt_MI,2),1)), squeeze(nanstd(nanmean(prefMaskOnly_noadapt_MI,2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefMaskOnly_singadapt_MI,2),1)), squeeze(nanstd(nanmean(prefMaskOnly_singadapt_MI,2),1)./sqrt(n_singadapt)), '-o')
title('All Mask')
xlabel('Test contrast')
ylabel('MI')
ylim([-1 1])
hline(0)
subplot(3,2,6)
n_singadapt = squeeze(sum(~isnan(nanmean(prefMaskOnly_singadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefMaskOnly_singadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefMaskOnly_singadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2),[],1))./sqrt(n_singadapt)),'-o');
title('All Mask')
xlabel('Test contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
hline(0)
legend({'No adapt', 'Single'})
suptitle(['Mask Pref- n = ' num2str(min(n_singadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_MI_prefMaskAll.pdf'),'-dpdf','-bestfit');

%Mask Resp
figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(respMask_noadapt_MI(:,it,:)),1);
    n_singadapt = sum(~isnan(respMask_singadapt_MI(:,it,:)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(respMask_noadapt_MI(:,it,:),1)), squeeze(nanstd(respMask_noadapt_MI(:,it,:),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(respMask_singadapt_MI(:,it,:),1)), squeeze(nanstd(respMask_singadapt_MI(:,it,:),1)./sqrt(n_singadapt)), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('MI')
    ylim([-1 1])
    hline(0)
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(respMask_noadapt_MI,2)),1);
n_singadapt = sum(~isnan(nanmean(respMask_singadapt_MI,2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(respMask_noadapt_MI,2),1)), squeeze(nanstd(nanmean(respMask_noadapt_MI,2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(respMask_singadapt_MI,2),1)), squeeze(nanstd(nanmean(respMask_singadapt_MI,2),1)./sqrt(n_singadapt)), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('MI')
hline(0)
ylim([-1 1])
subplot(3,2,6)
n_singadapt = squeeze(sum(~isnan(nanmean(respMask_singadapt_MI,2)-nanmean(respMask_noadapt_MI,2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(respMask_singadapt_MI,2)-nanmean(respMask_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(respMask_singadapt_MI,2)-nanmean(respMask_noadapt_MI,2),[],1))./sqrt(n_singadapt)),'-o');
title('All Test')
xlabel('Mask contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
hline(0)
legend({'No adapt', 'Single'})
suptitle(['Mask Resp- n = ' num2str(min(n_singadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_MI_maskRespAll.pdf'),'-dpdf','-bestfit');

figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefPlaidOnly_noadapt_MI(:,it,:)),1);
    n_singadapt = sum(~isnan(prefPlaidOnly_singadapt_MI(:,it,:)),1);
    n_contadapt = sum(~isnan(prefPlaidOnly_singadapt_MI(:,it,:)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefPlaidOnly_noadapt_MI(:,it,:),1)), squeeze(nanstd(prefPlaidOnly_noadapt_MI(:,it,:),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefPlaidOnly_singadapt_MI(:,it,:),1)), squeeze(nanstd(prefPlaidOnly_singadapt_MI(:,it,:),1)./sqrt(n_singadapt)), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('MI')
    ylim([-1 1])
    hline(0)
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefPlaidOnly_noadapt_MI,2)),1);
n_singadapt = sum(~isnan(nanmean(prefPlaidOnly_singadapt_MI,2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefPlaidOnly_noadapt_MI,2),1)), squeeze(nanstd(nanmean(prefPlaidOnly_noadapt_MI,2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefPlaidOnly_singadapt_MI,2),1)), squeeze(nanstd(nanmean(prefPlaidOnly_singadapt_MI,2),1)./sqrt(n_singadapt)), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('MI')
ylim([-1 1])
hline(0)
subplot(3,2,6)
n_singadapt = squeeze(sum(~isnan(nanmean(prefPlaidOnly_singadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefPlaidOnly_singadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefPlaidOnly_singadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2),[],1))./sqrt(n_singadapt)),'-o');
title('All Test')
xlabel('Mask contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
hline(0)
legend({'No adapt', 'Single'})
suptitle(['Plaid Pref- n = ' num2str(min(n_singadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_MI_prefPlaidOnly.pdf'),'-dpdf','-bestfit');

%%
prefTestOnly_noadapt_MIP = zeros(length(prefTestInd),nTest,nMask);
prefMaskOnly_noadapt_MIP = zeros(length(prefMaskInd),nTest,nMask);
prefPlaidOnly_noadapt_MIP = zeros(length(prefPlaidInd),nTest,nMask);
prefTestOnly_singadapt_MIP = zeros(length(prefTestInd),nTest,nMask);
prefMaskOnly_singadapt_MIP = zeros(length(prefMaskInd),nTest,nMask);
prefPlaidOnly_singadapt_MIP = zeros(length(prefPlaidInd),nTest,nMask);
respMask_noadapt_MIP = zeros(length(respmask_ind_all),nTest,nMask);
respMask_singadapt_MIP = zeros(length(respmask_ind_all),nTest,nMask);


for it = 1:nTest
    for im = 1:nMask
        prefTestOnly_noadapt_MIP(:,it,im) = squeeze((prefTestOnly_noadapt_respCells(:,im,it)-(prefTestOnly_noadapt_respCells(:,1,it)+prefTestOnly_noadapt_respCells(:,im,1)))./...
            (prefTestOnly_noadapt_respCells(:,im,it)+(prefTestOnly_noadapt_respCells(:,1,it)+prefTestOnly_noadapt_respCells(:,im,1))));
        prefMaskOnly_noadapt_MIP(:,it,im) = squeeze((prefMaskOnly_noadapt_respCells(:,im,it)-(prefMaskOnly_noadapt_respCells(:,1,it)+prefMaskOnly_noadapt_respCells(:,im,1)))./...
            (prefMaskOnly_noadapt_respCells(:,im,it)+(prefMaskOnly_noadapt_respCells(:,1,it)+prefMaskOnly_noadapt_respCells(:,im,1))));
        prefPlaidOnly_noadapt_MIP(:,it,im) = squeeze((prefPlaidOnly_noadapt_respCells(:,im,it)-(prefPlaidOnly_noadapt_respCells(:,1,it)+prefPlaidOnly_noadapt_respCells(:,im,1)))./...
            (prefPlaidOnly_noadapt_respCells(:,im,it)+(prefPlaidOnly_noadapt_respCells(:,1,it)+prefPlaidOnly_noadapt_respCells(:,im,1))));
        respMask_noadapt_MIP(:,it,im) = squeeze((respMask_noadapt_respCells(:,im,it)-(respMask_noadapt_respCells(:,1,it)+respMask_noadapt_respCells(:,im,1)))./...
            (respMask_noadapt_respCells(:,im,it)+(respMask_noadapt_respCells(:,1,it)+respMask_noadapt_respCells(:,im,1))));

        prefTestOnly_singadapt_MIP(:,it,im) = squeeze((prefTestOnly_singadapt_respCells(:,im,it)-(prefTestOnly_singadapt_respCells(:,1,it)+prefTestOnly_singadapt_respCells(:,im,1)))./...
            (prefTestOnly_singadapt_respCells(:,im,it)+(prefTestOnly_singadapt_respCells(:,1,it)+prefTestOnly_singadapt_respCells(:,im,1))));
        prefMaskOnly_singadapt_MIP(:,it,im) = squeeze((prefMaskOnly_singadapt_respCells(:,im,it)-(prefMaskOnly_singadapt_respCells(:,1,it)+prefMaskOnly_singadapt_respCells(:,im,1)))./...
            (prefMaskOnly_singadapt_respCells(:,im,it)+(prefMaskOnly_singadapt_respCells(:,1,it)+prefMaskOnly_singadapt_respCells(:,im,1))));
        prefPlaidOnly_singadapt_MIP(:,it,im) = squeeze((prefPlaidOnly_singadapt_respCells(:,im,it)-(prefPlaidOnly_singadapt_respCells(:,1,it)+prefPlaidOnly_singadapt_respCells(:,im,1)))./...
            (prefPlaidOnly_singadapt_respCells(:,im,it)+(prefPlaidOnly_singadapt_respCells(:,1,it)+prefPlaidOnly_singadapt_respCells(:,im,1))));
        respMask_singadapt_MIP(:,it,im) = squeeze((respMask_singadapt_respCells(:,im,it)-(respMask_singadapt_respCells(:,1,it)+respMask_singadapt_respCells(:,im,1)))./...
            (respMask_singadapt_respCells(:,im,it)+(respMask_singadapt_respCells(:,1,it)+respMask_singadapt_respCells(:,im,1))));
    end
end

figure; 
subplot(2,2,1)
imagesc(squeeze(nanmean(prefTestOnly_noadapt_MIP(:,2:5,2:5),1))); clim([-1 1]); colormap('redblue'); title('Test Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
colorbar
axis square
subplot(2,2,2)
imagesc(squeeze(nanmean(prefMaskOnly_noadapt_MIP(:,2:5,2:5),1))); clim([-1 1]); colormap('redblue'); title('Mask Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
axis square
colorbar
subplot(2,2,3)
imagesc(squeeze(nanmean(prefPlaidOnly_noadapt_MIP(:,2:5,2:5),1))); clim([-1 1]); colormap('redblue'); title('Plaid Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
colorbar
axis square
suptitle({'No adapt', '[M+T]-[M]-[T]/[M+T]+[M]+[T]'})
print(fullfile(fnout, 'adaptSummary_plaidModulationByPref-noAdapt.pdf'),'-dpdf','-bestfit');


figure; 
subplot(2,2,1)
imagesc(squeeze(nanmean(prefTestOnly_singadapt_MIP(:,2:5,2:5),1))); clim([-1 1]); colormap('redblue'); title('Test Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
colorbar
axis square
subplot(2,2,2)
imagesc(squeeze(nanmean(prefMaskOnly_singadapt_MIP(:,2:5,2:5),1))); clim([-1 1]); colormap('redblue'); title('Mask Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
axis square
colorbar
subplot(2,2,3)
imagesc(squeeze(nanmean(prefPlaidOnly_singadapt_MIP(:,2:5,2:5),1))); clim([-1 1]); colormap('redblue'); title('Plaid Pref')
ylabel('Test Contrast')
set(gca,'Xtick',[1:4])
set(gca,'Ytick',[1:4])
xticklabels(chop(maskCons(2:5),2))
xlabel('Mask Contrast')
yticklabels(chop(testCons(2:5),2))
colorbar
axis square
suptitle({'Sing adapt', '[M+T]-[M]-[T]/[M+T]+[M]+[T]'})
print(fullfile(fnout, 'adaptSummary_plaidModulationByPref-SingAdapt.pdf'),'-dpdf','-bestfit');

%Test pref
figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefTestOnly_noadapt_MIP(:,it,2:5)),1);
    n_singadapt = sum(~isnan(prefTestOnly_singadapt_MIP(:,it,2:5)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOnly_noadapt_MIP(:,it,2:5),1)), squeeze(nanstd(prefTestOnly_noadapt_MIP(:,it,2:5),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOnly_singadapt_MIP(:,it,2:5),1)), squeeze(nanstd(prefTestOnly_singadapt_MIP(:,it,2:5),1)./sqrt(n_singadapt)), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('SI')
    ylim([-1 1])
    hline(0)
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefTestOnly_noadapt_MIP(:,2:5,2:5),2)),1);
n_singadapt = sum(~isnan(nanmean(prefTestOnly_singadapt_MIP(:,2:5,2:5),2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOnly_noadapt_MIP(:,2:5,2:5),2),1)), squeeze(nanstd(nanmean(prefTestOnly_noadapt_MIP(:,2:5,2:5),2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOnly_singadapt_MIP(:,2:5,2:5),2),1)), squeeze(nanstd(nanmean(prefTestOnly_singadapt_MIP(:,2:5,2:5),2),1)./sqrt(n_singadapt)), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('SI')
ylim([-1 1])
hline(0)
subplot(3,2,6)
n_singadapt = squeeze(sum(~isnan(nanmean(prefTestOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefTestOnly_noadapt_MIP(:,2:5,2:5),2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefTestOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefTestOnly_noadapt_MIP(:,2:5,2:5),2),1)), (squeeze(nanstd(nanmean(prefTestOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefTestOnly_noadapt_MIP(:,2:5,2:5),2),[],1))./sqrt(n_singadapt)),'-o');
title('All Test')
xlabel('Mask contrast')
ylabel('MIP (post-pre)')
ylim([-1 1])
hline(0)
legend({'No adapt', 'Single'})
suptitle(['Test Pref- n = ' num2str(min(n_singadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_SI_prefTestAll.pdf'),'-dpdf','-bestfit');
%Mask pref
figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefMaskOnly_noadapt_MIP(:,it,2:5)),1);
    n_singadapt = sum(~isnan(prefMaskOnly_singadapt_MIP(:,it,2:5)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefMaskOnly_noadapt_MIP(:,it,2:5),1)), squeeze(nanstd(prefMaskOnly_noadapt_MIP(:,it,2:5),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefMaskOnly_singadapt_MIP(:,it,2:5),1)), squeeze(nanstd(prefMaskOnly_singadapt_MIP(:,it,2:5),1)./sqrt(n_singadapt)), '-o')
    title(['Mask = ' num2str(testCons(it))])
    xlabel('Test contrast')
    ylabel('SI')
    ylim([-1 1])
    hline(0)
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefMaskOnly_noadapt_MIP(:,2:5,2:5),2)),1);
n_singadapt = sum(~isnan(nanmean(prefMaskOnly_singadapt_MIP(:,2:5,2:5),2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefMaskOnly_noadapt_MIP(:,2:5,2:5),2),1)), squeeze(nanstd(nanmean(prefMaskOnly_noadapt_MIP(:,2:5,2:5),2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefMaskOnly_singadapt_MIP(:,2:5,2:5),2),1)), squeeze(nanstd(nanmean(prefMaskOnly_singadapt_MIP(:,2:5,2:5),2),1)./sqrt(n_singadapt)), '-o')
title('All Mask')
xlabel('Test contrast')
ylabel('SI')
ylim([-1 1])
hline(0)
subplot(3,2,6)
n_singadapt = squeeze(sum(~isnan(nanmean(prefMaskOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefMaskOnly_noadapt_MIP(:,2:5,2:5),2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefMaskOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefMaskOnly_noadapt_MIP(:,2:5,2:5),2),1)), (squeeze(nanstd(nanmean(prefMaskOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefMaskOnly_noadapt_MIP(:,2:5,2:5),2),[],1))./sqrt(n_singadapt)),'-o');
title('All Mask')
xlabel('Test contrast')
ylabel('MIP (post-pre)')
ylim([-1 1])
hline(0)
legend({'No adapt', 'Single'})
suptitle(['Mask Pref- n = ' num2str(min(n_singadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_SI_prefMaskAll.pdf'),'-dpdf','-bestfit');

figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefPlaidOnly_noadapt_MIP(:,it,2:5)),1);
    n_singadapt = sum(~isnan(prefPlaidOnly_singadapt_MIP(:,it,2:5)),1);
    n_contadapt = sum(~isnan(prefPlaidOnly_singadapt_MIP(:,it,2:5)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefPlaidOnly_noadapt_MIP(:,it,2:5),1)), squeeze(nanstd(prefPlaidOnly_noadapt_MIP(:,it,2:5),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefPlaidOnly_singadapt_MIP(:,it,2:5),1)), squeeze(nanstd(prefPlaidOnly_singadapt_MIP(:,it,2:5),1)./sqrt(n_singadapt)), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('SI')
    ylim([-1 1])
    hline(0)
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefPlaidOnly_noadapt_MIP(:,2:5,2:5),2)),1);
n_singadapt = sum(~isnan(nanmean(prefPlaidOnly_singadapt_MIP(:,2:5,2:5),2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefPlaidOnly_noadapt_MIP(:,2:5,2:5),2),1)), squeeze(nanstd(nanmean(prefPlaidOnly_noadapt_MIP(:,2:5,2:5),2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefPlaidOnly_singadapt_MIP(:,2:5,2:5),2),1)), squeeze(nanstd(nanmean(prefPlaidOnly_singadapt_MIP(:,2:5,2:5),2),1)./sqrt(n_singadapt)), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('SI')
ylim([-1 1])
hline(0)
subplot(3,2,6)
n_singadapt = squeeze(sum(~isnan(nanmean(prefPlaidOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefPlaidOnly_noadapt_MIP(:,2:5,2:5),2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefPlaidOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefPlaidOnly_noadapt_MIP(:,2:5,2:5),2),1)), (squeeze(nanstd(nanmean(prefPlaidOnly_singadapt_MIP(:,2:5,2:5),2)-nanmean(prefPlaidOnly_noadapt_MIP(:,2:5,2:5),2),[],1))./sqrt(n_singadapt)),'-o');
title('All Test')
xlabel('Mask contrast')
ylabel('MIP (post-pre)')
ylim([-1 1])
hline(0)
legend({'No adapt', 'Single'})
suptitle(['Plaid Pref- n = ' num2str(min(n_singadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_SI_prefPlaidOnly.pdf'),'-dpdf','-bestfit');

%% histograms
figure;
for im = 2:nMask
    subplot(2,2,im-1)
    hist(respMask_noadapt_MIP(:,im,im)-respMask_singadapt_MIP(:,im,im))
    xlabel('NoAdapt-Adapt MIP')
    title(['Test/Mask =  ' num2str(maskCons(im))])
    vline(0)
end

figure;
for im = 2:nMask
    subplot(2,2,im-1)
    hist(prefMaskOnly_noadapt_MIP(:,im,im)-prefMaskOnly_singadapt_MIP(:,im,im))
    xlabel('NoAdapt-Adapt MIP')
    title(['Test/Mask =  ' num2str(maskCons(im))])
    vline(nanmean(prefMaskOnly_noadapt_MIP(:,im,im)-prefMaskOnly_singadapt_MIP(:,im,im),1))
end


figure;
for im = 2:nMask
    subplot(2,4,im-1)
    hist(prefMaskOnly_noadapt_MIP(:,im,im))
    xlabel('No Adapt MIP')
    title(['Test/Mask =  ' num2str(maskCons(im))])
    vline(nanmean(prefMaskOnly_noadapt_MIP(:,im,im),1))
    xlim([-1 1])
    subplot(2,4,im-1+4)
    hist(prefMaskOnly_singadapt_MIP(:,im,im))
    xlabel('Adapt MIP')
    title(['Test/Mask =  ' num2str(maskCons(im))])
    vline(nanmean(prefMaskOnly_singadapt_MIP(:,im,im),1))
    xlim([-1 1])
end

figure;
for im = 2:nMask
    subplot(2,2,im-1)
    cdfplot(prefMaskOnly_noadapt_MIP(:,im,im))
    hold on
    cdfplot(prefMaskOnly_singadapt_MIP(:,im,im))
    xlabel('Suppression Index')
    title(['Test/Mask =  ' num2str(maskCons(im))])
    xlim([-1 1])
end
legend({'No adapt', 'Single'}, 'location','southeast')
suptitle('Mask pref')
print(fullfile(fnout, 'adaptSummary_SIcdfs_prefMaskOnly.pdf'),'-dpdf','-bestfit');

figure;
for im = 2:nMask
    subplot(2,2,im-1)
    cdfplot(prefPlaidOnly_noadapt_MIP(:,im,im))
    hold on
    cdfplot(prefPlaidOnly_singadapt_MIP(:,im,im))
    xlabel('No Adapt MIP')
    title(['Test/Mask =  ' num2str(maskCons(im))])
    xlim([-1 1])
end
legend({'No adapt', 'Single'}, 'location','northwest')
suptitle('Plaid pref')

figure;
for im = 2:nMask
    subplot(2,2,im-1)
    cdfplot(respMask_noadapt_MIP(:,im,im))
    hold on
    cdfplot(respMask_singadapt_MIP(:,im,im))
    xlabel('No Adapt MIP')
    title(['Test/Mask =  ' num2str(maskCons(im))])
    xlim([-1 1])
end
legend({'No adapt', 'Single'}, 'location','northwest')
suptitle('Mask resp')

figure;
ind_l = find(respMask_noadapt_MIP(:,5,5)<-0.2);
nC_l = length(ind_l);
ind_h = find(respMask_noadapt_MIP(:,5,5)>0.2);
nC_h = length(ind_h);
plot(repmat([1 2], [nC_l 1])',[respMask_noadapt_MIP(ind_l,5,5) respMask_singadapt_MIP(ind_l,5,5)]','b')
hold on
errorbar([1; 2], [nanmean(respMask_noadapt_MIP(ind_l,5,5),1) nanmean(respMask_singadapt_MIP(ind_l,5,5),1)]', [nanstd(respMask_noadapt_MIP(ind_l,5,5),[],1)./sqrt(nC_l) nanstd(respMask_singadapt_MIP(ind_l,5,5),[],1)./sqrt(nC_l)]','ob')
plot(repmat([1 2], [nC_h 1])',[respMask_noadapt_MIP(ind_h,5,5) respMask_singadapt_MIP(ind_h,5,5)]','r')
errorbar([1; 2], [nanmean(respMask_noadapt_MIP(ind_h,5,5),1) nanmean(respMask_singadapt_MIP(ind_h,5,5),1)]', [nanstd(respMask_noadapt_MIP(ind_h,5,5),[],1)./sqrt(nC_h) nanstd(respMask_singadapt_MIP(ind_h,5,5),[],1)./sqrt(nC_h)]','or')
xlim([0 3])

figure;
respMask_MIP_diff = respMask_noadapt_MIP-respMask_singadapt_MIP;
cdfplot(respMask_MIP_diff(ind_l,5,5))
hold on
cdfplot(respMask_MIP_diff(ind_h,5,5))