%% current datasets
mouse_mat = strvcat('i674', 'i674', 'i684', 'i696');
date_mat = strvcat('170112', '170116','170207','170210');
run_str_mat = strvcat('runs-002-003','runs-003-004','runs-004-005','runs-002-003');

nexp = size(mouse_mat,1);
mouse_str = [];
for iexp = 1:nexp
    mouse_str = [mouse_str mouse_mat(iexp,:)];
end

%% rerun datasets
nexp = size(mouse_mat,1);
ind_min = cell(1,nexp);
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    run_str = run_str_mat(iexp,:);
    load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))
    load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfofData.mat']))
    FS_RandInt_ROC
end
%% collect datasets
nexp = size(mouse_mat,1);
good_ind_base_all = [];
good_ind_targ_all = [];
roc_base_all = [];
roc_base_N_all = [];
roc_base_N1_all = [];
roc_targ_all = [];
roc_targ_N_all = [];
roc_targ_N1_all = [];
targ_resp_N1_all = [];
base_resp_N1_all = [];
targ_resp_all = [];
base_resp_all = [];
baseresp_norm_all = [];
baseresp_ioff_norm_all = [];
baseresp_ibin_norm_all = [];
nCells = 0;
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    run_str = run_str_mat(iexp,:);
    load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_ROC.mat']))
    roc_base_all = [roc_base_all roc_base];
    roc_base_N_all = [roc_base_N_all roc_base_allN];
    roc_base_N1_all = cat(3,roc_base_N1_all, roc_base_N1);
    roc_targ_all = cat(3,roc_targ_all, roc_targ);
    roc_targ_N_all = cat(3,roc_targ_N_all, roc_targ_allN);
    roc_targ_N1_all = cat(4,roc_targ_N1_all, roc_targ_N1);
    base_resp_all = cat(3,base_resp_all, base_resp);
    targ_resp_all = cat(4,targ_resp_all, targ_resp);
    base_resp_N1_all = cat(3,base_resp_N1_all, base_resp_N1);
    targ_resp_N1_all = cat(4,targ_resp_N1_all, targ_resp_N1);
    good_ind_base_all = [good_ind_base_all good_ind_base+nCells];
    good_ind_targ_all = [good_ind_targ_all good_ind_targ+nCells];
    nCells = size(roc_base_all, 2);
    load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_baseResp.mat']))
    baseresp_norm_all = [baseresp_norm_all; baseresp_dfof_norm];
    baseresp_ioff_norm_all = cat(1,baseresp_ioff_norm_all, baseresp_dfof_ioff_norm);
    baseresp_ibin_norm_all = [baseresp_ibin_norm_all; baseresp_dfof_ibin_norm];
end
load(fullfile(LG_base, '\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_stimData.mat']))

%% plot figures
ISIs = [250 500 750];
figure; 
subplot(1,3,1)
errorbar(ISIs,mean(roc_base_all(:,good_ind_base_all),2),std(roc_base_all(:,good_ind_base_all),[],2)./sqrt(length(good_ind_base_all)),'-ok')
axis square
ylim([0.3 0.7])
xlabel('ISI (ms)')
ylabel('auROC')
hline(0.5)
title('Base- All N-1')
for i = 1:nDelta
    subplot(1,3,i+1)
    errorbar(ISIs,mean(roc_targ_all(:,i,good_ind_targ_all),3),std(roc_targ_all(:,i,good_ind_targ_all),[],3)./sqrt(length(good_ind_targ_all)),'-ok')
    axis square
    ylim([0.4 0.6])
    xlabel('ISI (ms)')
    ylabel('auROC')
    hline(0.5)
    title([num2str(deltas(i)) ' deg - All N-1'])
end
suptitle([mouse_str '- All N-1 auROC'])
print(fullfile(LG_base, '\Analysis\2P', 'Adaptation',['randInt_allNminus1_summary.pdf']),'-dpdf','-fillpage')

figure; 
subplot(1,3,1)
errorbar(ISIs,mean(roc_base_N_all(:,good_ind_base_all),2),std(roc_base_N_all(:,good_ind_base_all),[],2)./sqrt(length(good_ind_base_all)),'-ok')
axis square
ylim([0.3 0.7])
xlabel('ISI (ms)')
ylabel('auROC')
hline(0.5)
title('Base- All N')
for i = 1:nDelta
    subplot(1,3,i+1)
    errorbar(ISIs,mean(roc_targ_N_all(:,i,good_ind_targ_all),3),std(roc_targ_N_all(:,i,good_ind_targ_all),[],3)./sqrt(length(good_ind_targ_all)),'-ok')
    axis square
    ylim([0.4 0.6])
    xlabel('ISI (ms)')
    ylabel('auROC')
    hline(0.5)
    title([num2str(deltas(i)) ' deg - All N'])
end
suptitle([mouse_str '- All N auROC'])
print(fullfile(LG_base, '\Analysis\2P', 'Adaptation',['randInt_allN_summary.pdf']),'-dpdf','-fillpage')
[p_roc_base table_roc_base stats_roc_base] = anova1(roc_base_N_all(:,good_ind_base_all)');
[p_roc_targ30 table_roc_targ30 stats_roc_targ30] = anova1(squeeze(roc_targ_N_all(:,1,good_ind_base_all))');

figure;
roc_base_N1_avg = zeros(noff,noff,2);
p_roc_base_N1 = zeros(noff,2);
start = 1;
for i = 1:noff
    subplot(3,3,start)
    errorbar(ISIs,squeeze(mean(roc_base_N1_all(i,:,good_ind_base_all),3)),squeeze(std(roc_base_N1_all(i,:,good_ind_base_all),[],3))./sqrt(length(good_ind_base_all)),'-ok')
    roc_base_N1_avg(i,:,1) = mean(roc_base_N1_all(i,:,good_ind_base_all),3);
    roc_base_N1_avg(i,:,2) = std(roc_base_N1_all(i,:,good_ind_base_all),[],3)./sqrt(length(good_ind_base_all));
    [h p_roc_base_N1(i,1)] = ttest(roc_base_N1_all(i,2,good_ind_base_all),roc_base_N1_all(i,1,good_ind_base_all));
    [h p_roc_base_N1(i,2)] = ttest(roc_base_N1_all(i,3,good_ind_base_all),roc_base_N1_all(i,1,good_ind_base_all));
    axis square
    ylim([0.4 0.6])
    xlabel('ISI (ms)')
    ylabel('auROC')
    hline(0.5)
    title(['Base - N = ' num2str(chop(offs(i)*(1000/frameRateHz),3)) ' ms'])
    for idelta = 1:nDelta
        subplot(3,3,start+idelta)
        errorbar(ISIs,squeeze(mean(roc_targ_N1_all(i,:,idelta,good_ind_targ_all),4)),squeeze(std(roc_targ_N1_all(i,:,idelta,good_ind_targ_all),[],4))./sqrt(length(good_ind_targ_all)),'-ok')
        axis square
        ylim([0.4 0.6])
        xlabel('ISI (ms)')
        ylabel('auROC')
        hline(0.5)
        title([num2str(deltas(idelta)) ' deg - N = ' num2str(chop(offs(i)*(1000/frameRateHz),3)) ' ms'])
    end
    start = start+3;
end
suptitle([mouse_str '- auROC by N-1'])
print(fullfile(LG_base, '\Analysis\2P', 'Adaptation',['randInt_byNminus1_summary.pdf']),'-dpdf','-fillpage')

figure;
start = 1;
for i = 1:noff
    subplot(4,2,start)
    errorbar(ISIs, squeeze(mean(base_resp_all(i,:,good_ind_base_all),3)), squeeze(std(base_resp_all(i,:,good_ind_base_all),[],3)./sqrt(length(good_ind_base_all))), '-ok')
    title(['N = ' num2str(chop(offs(i).*(1000/frameRateHz),3)) '; N response'])
    xlabel('N-1 ISI')
    ylabel('Normalized dF/F')
    ylim([0 0.8])
    subplot(4,2,start+1)
    errorbar(ISIs, squeeze(mean(base_resp_N1_all(i,:,good_ind_base_all),3)), squeeze(std(base_resp_N1_all(i,:,good_ind_base_all),[],3)./sqrt(length(good_ind_base_all))), '-ok')
    start = start+2;
    title(['N = ' num2str(chop(offs(i).*(1000/frameRateHz),3)) '; N-1 response'])
    xlabel('N-1 ISI')
    ylabel('Normalized dF/F')
    ylim([0 0.8])
end
subplot(4,2,7)
errorbar(ISIs, squeeze(mean(mean(base_resp_all(:,:,good_ind_base_all),1),3)), squeeze(std(mean(base_resp_all(:,:,good_ind_base_all),1),[],3)./sqrt(length(good_ind_base_all))), '-ok')
base_resp_all_avg(:,:,1) = squeeze(mean(mean(base_resp_all(:,:,good_ind_base_all),1),3));
base_resp_all_avg(:,:,2) = squeeze(std(mean(base_resp_all(:,:,good_ind_base_all),1),[],3)./sqrt(length(good_ind_base_all)));
title(['All N- N response'])
xlabel('N-1 ISI')
ylabel('Normalized dF/F')
ylim([0 0.8])
subplot(4,2,8)
errorbar(ISIs, squeeze(mean(mean(base_resp_N1_all(:,:,good_ind_base_all),1),3)), squeeze(std(mean(base_resp_N1_all(:,:,good_ind_base_all),1),[],3)./sqrt(length(good_ind_base_all))), '-ok')
base_resp_N1_all_avg(:,:,1) = squeeze(mean(mean(base_resp_N1_all(:,:,good_ind_base_all),1),3));
base_resp_N1_all_avg(:,:,2) = squeeze(std(mean(base_resp_N1_all(:,:,good_ind_base_all),1),[],3)./sqrt(length(good_ind_base_all)));
title(['All N- N-1 response'])
xlabel('N-1 ISI')
ylabel('Normalized dF/F')
ylim([0 0.8])
suptitle([mouse_str '- Norm dF/F by N-1'])
print(fullfile(LG_base, '\Analysis\2P', 'Adaptation',['randInt_dFoverF_byNminus1Int.pdf']),'-dpdf','-fillpage')

p_allN_Nresp = anova1(squeeze(mean(base_resp_all(:,:,good_ind_base_all),1))');
p_allN_N1resp = anova1(squeeze(mean(base_resp_N1_all(:,:,good_ind_base_all),1))');
%% Trial length analysis
figure;
%average PP ratio by cycle
subplot(2,2,1)
errorbar(1:5, mean(baseresp_norm_all(good_ind_base_all,:),1), std(baseresp_norm_all(good_ind_base_all,:),[],1)./sqrt(length(good_ind_base_all)),'-o')
xlabel('Cycle #')
ylabel('Normalized amplitude')
ylim([0 1])
title(['n = ' num2str(length(good_ind_base_all))])
%average PP ratio by length
subplot(2,2,2)
n = sum(~isnan(baseresp_ibin_norm_all(good_ind_base_all,:)),1);
errorbar(x, nanmean(baseresp_ibin_norm_all(good_ind_base_all,:),1), nanstd(baseresp_ibin_norm_all(good_ind_base_all,:),[],1)./sqrt(n),'-o')
xlabel('Trial length (ms)')
ylabel('Normalized amplitude')
ylim([0 1])
xlim([-500 4000])
%average PP ratio by interval
subplot(2,2,3)
for i = 1:noff
    n = sum(~isnan(baseresp_ioff_norm_all(good_ind_base_all,:,i)),1);
    errorbar(1:5, nanmean(baseresp_ioff_norm_all(good_ind_base_all,:,i),1), nanstd(baseresp_ioff_norm_all(good_ind_base_all,:,i),[],1)./sqrt(n),'-o')
    hold on
end
xlabel('Cycle #')
ylabel('Normalized amplitude')
ylim([0 1])
title('By interval')
suptitle([mouse_str '- Norm dF/F trial length'])
print(fullfile(LG_base, '\Analysis\2P', 'Adaptation',['randInt_dFoverF_trialLength.pdf']),'-dpdf','-fillpage')

%%
save(fullfile(LG_base, '\Analysis\2P', 'Adaptation', 'randIntSummary.mat'),'base_resp_N1_all', 'targ_resp_N1_all', 'roc_base_all', 'roc_base_N1_all', 'roc_targ_all', 'roc_targ_N1_all', 'good_ind_base_all', 'good_ind_targ_all')