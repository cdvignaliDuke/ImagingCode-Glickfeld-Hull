%contingent adaptation experiments
%250 ms adapt stim
mouse_mat = strvcat('i1303','i1103','i1304','i1303','i1304','i1103','i1303');
date_mat = strvcat('190501','190508','190508','190619','190619','190621','190621');
run_mat = strvcat('002','002','002','003','003','002','002');

%1s adapt stim
% mouse_mat = strvcat('i1303','i1304','i1103','i1303');
% date_mat = strvcat('190810','190817','190824','190824');
% run_mat = strvcat('002','002','002','002');
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';
nexp = size(mouse_mat,1);

%% collect datasets
totCells = 0;
preftestonly_ind_all = [];
prefmaskonly_ind_all = [];
prefplaidonly_ind_all = [];
resptest_ind_all = [];
respmask_ind_all = [];
respplaid_ind_all = [];
noadapt_resp_cell_all= cell(5,5);
asynadapt_resp_cell_all= cell(5,5);
contadapt_resp_cell_all= cell(5,5);
for iexp = 1:nexp
    mouse = mouse_mat(iexp,:);
    date = date_mat(iexp,:);
    run_str = catRunName(run_mat(iexp,:), 1);
    
    fn = fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]);
    load(fullfile(fn,[date '_' mouse '_' run_str '_allCellResp.mat']));
    %load(fullfile(fn,[date '_' mouse '_' run_str '_allCellRespLate.mat']));
    load(fullfile(fn,[date '_' mouse '_' run_str '_dataStim.mat']));
    nMask = length(maskCons);
    nTest = length(testCons);

    for im = 1:nMask
        for it = 1:nTest
            noadapt_resp_cell_all{im,it} = [noadapt_resp_cell_all{im,it}; nanmean(noadapt_resp_cell{im,it},2)];
            asynadapt_resp_cell_all{im,it} = [asynadapt_resp_cell_all{im,it}; nanmean(asynadapt_resp_cell{im,it},2)];
            contadapt_resp_cell_all{im,it} = [contadapt_resp_cell_all{im,it}; nanmean(contadapt_resp_cell{im,it},2)];
        end
    end
    preftestonly_ind_all = [preftestonly_ind_all; preftestonly_ind+totCells];
    prefmaskonly_ind_all = [prefmaskonly_ind_all; prefmaskonly_ind+totCells];
    prefplaidonly_ind_all = [prefplaidonly_ind_all; prefplaidonly_ind+totCells];
    resptest_ind_all = [resptest_ind_all; resptest_ind+totCells];
    respmask_ind_all = [respmask_ind_all; respmask_ind+totCells];
    respplaid_ind_all = [respplaid_ind_all; respplaid_ind+totCells];
    
    totCells = totCells+size(noadapt_resp_cell{1,1},1);
    fprintf([num2str(length(preftestonly_ind)) ' ' num2str(length(prefmaskonly_ind)) ' ' num2str(length(prefplaidonly_ind)) '\n'])
end

%% 
fnout = fullfile(LG_base, 'Analysis\2P\ContingentAdaptation');
if ~exist(fnout)
    mkdir(fnout);
end

prefTestOnly_noadapt_resp = zeros(nMask,nTest,2);
prefMaskOnly_noadapt_resp = zeros(nMask,nTest,2);
prefPlaidOnly_noadapt_resp = zeros(nMask,nTest,2);
prefTestOnly_contadapt_resp = zeros(nMask,nTest,2);
prefMaskOnly_contadapt_resp = zeros(nMask,nTest,2);
prefPlaidOnly_contadapt_resp = zeros(nMask,nTest,2);
prefTestOnly_asynadapt_resp = zeros(nMask,nTest,2);
prefMaskOnly_asynadapt_resp = zeros(nMask,nTest,2);
prefPlaidOnly_asynadapt_resp = zeros(nMask,nTest,2);
prefTestOrMask_noadapt_resp = zeros(nMask,nTest,2);
prefTestOrMask_contadapt_resp = zeros(nMask,nTest,2);
prefTestOrMask_asynadapt_resp = zeros(nMask,nTest,2);

prefTestInd = intersect(preftestonly_ind_all,resptest_ind_all);
prefMaskInd = intersect(prefmaskonly_ind_all,respmask_ind_all);
%prefTestInd = intersect(find(noadapt_resp_cell_all{1,5}>noadapt_resp_cell_all{5,1}),resptest_ind_all);
%prefMaskInd = intersect(find(noadapt_resp_cell_all{1,5}<noadapt_resp_cell_all{5,1}),respmask_ind_all);
prefPlaidInd = intersect(prefplaidonly_ind_all,respplaid_ind_all);

for im = 1:nMask
    for it = 1:nTest
        prefTestOnly_noadapt_resp(im,it,1) = nanmean(noadapt_resp_cell_all{im,it}(prefTestInd,:),1);
        prefTestOnly_noadapt_resp(im,it,2) = nanstd(noadapt_resp_cell_all{im,it}(prefTestInd,:),[],1)./sqrt(length(prefTestInd));
        prefMaskOnly_noadapt_resp(im,it,1) = nanmean(noadapt_resp_cell_all{im,it}(prefMaskInd,:),1);
        prefMaskOnly_noadapt_resp(im,it,2) = nanstd(noadapt_resp_cell_all{im,it}(prefMaskInd,:),[],1)./sqrt(length(prefMaskInd));
        prefPlaidOnly_noadapt_resp(im,it,1) = nanmean(noadapt_resp_cell_all{im,it}(prefPlaidInd,:),1);
        prefPlaidOnly_noadapt_resp(im,it,2) = nanstd(noadapt_resp_cell_all{im,it}(prefPlaidInd,:),[],1)./sqrt(length(prefPlaidInd));
        prefTestOnly_contadapt_resp(im,it,1) = nanmean(contadapt_resp_cell_all{im,it}(prefTestInd,:),1);
        prefTestOnly_contadapt_resp(im,it,2) = nanstd(contadapt_resp_cell_all{im,it}(prefTestInd,:),[],1)./sqrt(length(prefTestInd));
        prefMaskOnly_contadapt_resp(im,it,1) = nanmean(contadapt_resp_cell_all{im,it}(prefMaskInd,:),1);
        prefMaskOnly_contadapt_resp(im,it,2) = nanstd(contadapt_resp_cell_all{im,it}(prefMaskInd,:),[],1)./sqrt(length(prefMaskInd));
        prefPlaidOnly_contadapt_resp(im,it,1) = nanmean(contadapt_resp_cell_all{im,it}(prefPlaidInd,:),1);
        prefPlaidOnly_contadapt_resp(im,it,2) = nanstd(contadapt_resp_cell_all{im,it}(prefPlaidInd,:),[],1)./sqrt(length(prefPlaidInd));
        prefTestOnly_asynadapt_resp(im,it,1) = nanmean(asynadapt_resp_cell_all{im,it}(prefTestInd,:),1);
        prefTestOnly_asynadapt_resp(im,it,2) = nanstd(asynadapt_resp_cell_all{im,it}(prefTestInd,:),[],1)./sqrt(length(prefTestInd));
        prefMaskOnly_asynadapt_resp(im,it,1) = nanmean(asynadapt_resp_cell_all{im,it}(prefMaskInd,:),1);
        prefMaskOnly_asynadapt_resp(im,it,2) = nanstd(asynadapt_resp_cell_all{im,it}(prefMaskInd,:),[],1)./sqrt(length(prefMaskInd));
        prefPlaidOnly_asynadapt_resp(im,it,1) = nanmean(asynadapt_resp_cell_all{im,it}(prefPlaidInd,:),1);
        prefPlaidOnly_asynadapt_resp(im,it,2) = nanstd(asynadapt_resp_cell_all{im,it}(prefPlaidInd,:),[],1)./sqrt(length(prefPlaidInd));
        prefTestOrMask_noadapt_temp = [noadapt_resp_cell_all{im,it}(prefTestInd,:); noadapt_resp_cell_all{it,im}(prefMaskInd,:)];
        prefTestOrMask_contadapt_temp = [contadapt_resp_cell_all{im,it}(prefTestInd,:); contadapt_resp_cell_all{it,im}(prefMaskInd,:)];
        prefTestOrMask_asynadapt_temp = [asynadapt_resp_cell_all{im,it}(prefTestInd,:); asynadapt_resp_cell_all{it,im}(prefMaskInd,:)];
        prefTestOrMask_noadapt_resp(im,it,1) = nanmean(prefTestOrMask_noadapt_temp,1);
        prefTestOrMask_noadapt_resp(im,it,2) = nanstd(prefTestOrMask_noadapt_temp,[],1)./sqrt(size(prefTestOrMask_noadapt_temp,1));
        prefTestOrMask_contadapt_resp(im,it,1) = nanmean(prefTestOrMask_contadapt_temp,1);
        prefTestOrMask_contadapt_resp(im,it,2) = nanstd(prefTestOrMask_contadapt_temp,[],1)./sqrt(size(prefTestOrMask_noadapt_temp,1));
        prefTestOrMask_asynadapt_resp(im,it,1) = nanmean(prefTestOrMask_asynadapt_temp,1);
        prefTestOrMask_asynadapt_resp(im,it,2) = nanstd(prefTestOrMask_asynadapt_temp,[],1)./sqrt(size(prefTestOrMask_asynadapt_temp,1));
    end
end

% TestOnly
figure;
subplot(1,3,1)
for im = 1:nMask
    errorbar(testCons, prefTestOnly_noadapt_resp(im,:,1),prefTestOnly_noadapt_resp(im,:,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,2)
for im = 1:nMask
    errorbar(testCons, prefTestOnly_contadapt_resp(im,:,1),prefTestOnly_contadapt_resp(im,:,2),'-o')
    hold on
end
title('Contingent')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,3)
for im = 1:nMask
    errorbar(testCons, prefTestOnly_asynadapt_resp(im,:,1),prefTestOnly_asynadapt_resp(im,:,2),'-o')
    hold on
end
title('Asynchronous')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend(num2str(maskCons'))
suptitle(['Test preferring cells- n = ' num2str(length(prefTestInd))])
print(fullfile(fnout, 'adaptSummary_testPrefAll_250ms.pdf'),'-dpdf','-bestfit');

%Mask Only
figure;
subplot(1,3,1)
for it = 1:nTest
    errorbar(maskCons, prefMaskOnly_noadapt_resp(:,it,1),prefMaskOnly_noadapt_resp(:,it,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,2)
for it = 1:nTest
    errorbar(maskCons, prefMaskOnly_contadapt_resp(:,it,1),prefMaskOnly_contadapt_resp(:,it,2),'-o')
    hold on
end
title('Contingent')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,3)
for it = 1:nTest
    errorbar(maskCons, prefMaskOnly_asynadapt_resp(:,it,1),prefMaskOnly_asynadapt_resp(:,it,2),'-o')
    hold on
end
title('Asynchronous')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend(num2str(testCons'))
suptitle(['Mask preferring cells- n = ' num2str(length(prefMaskInd))])
print(fullfile(fnout, 'adaptSummary_maskPrefAll_250ms.pdf'),'-dpdf','-bestfit');

% Plaid Only
figure;
subplot(1,3,1)
for im = 1:nMask
    errorbar(testCons, prefPlaidOnly_noadapt_resp(im,:,1),prefPlaidOnly_noadapt_resp(im,:,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,2)
for im = 1:nMask
    errorbar(testCons, prefPlaidOnly_contadapt_resp(im,:,1),prefPlaidOnly_contadapt_resp(im,:,2),'-o')
    hold on
end
title('Contingent')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,3)
for im = 1:nMask
    errorbar(testCons, prefPlaidOnly_asynadapt_resp(im,:,1),prefPlaidOnly_asynadapt_resp(im,:,2),'-o')
    hold on
end
title('Asynchronous')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend(num2str(maskCons'))
suptitle(['Plaid preferring cells- n = ' num2str(length(prefPlaidInd))])
print(fullfile(fnout, 'adaptSummary_plaidPrefOnly_250ms.pdf'),'-dpdf','-bestfit');

%Test or Mask
figure;
subplot(1,3,1)
for im = 1:nMask
    errorbar(testCons, prefTestOrMask_noadapt_resp(im,:,1),prefTestOrMask_noadapt_resp(im,:,2),'-o')
    hold on
end
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,2)
for im = 1:nMask
    errorbar(testCons, prefTestOrMask_contadapt_resp(im,:,1),prefTestOrMask_contadapt_resp(im,:,2),'-o')
    hold on
end
title('Contingent')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,3)
for im = 1:nMask
    errorbar(testCons, prefTestOrMask_asynadapt_resp(im,:,1),prefTestOrMask_asynadapt_resp(im,:,2),'-o')
    hold on
end
title('Asynchronous')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend(num2str(maskCons'))
suptitle(['Test or Mask preferring cells- n = ' num2str(length([prefTestInd; prefMaskInd]))])
print(fullfile(fnout, 'adaptSummary_prefTestOrMaskAll_250ms.pdf'),'-dpdf','-bestfit');

%by mask
figure;
for im = 1:nMask
    subplot(2,3,im)
    errorbar(testCons, prefTestOrMask_noadapt_resp(im,:,1),prefTestOrMask_noadapt_resp(im,:,2),'-o')
    hold on
    errorbar(testCons, prefTestOrMask_contadapt_resp(im,:,1),prefTestOrMask_contadapt_resp(im,:,2),'-o')
    errorbar(testCons, prefTestOrMask_asynadapt_resp(im,:,1),prefTestOrMask_asynadapt_resp(im,:,2),'-o')
    title(['Mask = ' num2str(maskCons(im))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
end
legend('No adapt', 'Contingent', 'Asynchronous');
suptitle(['Test or Mask preferring cells- n = ' num2str(length([prefTestInd; prefMaskInd]))])
print(fullfile(fnout, 'adaptSummary_byMask_prefTestOrMaskAll_250ms.pdf'),'-dpdf','-bestfit');

figure;
for im = 1:nMask
    subplot(2,3,im)
    errorbar(testCons, prefTestOnly_noadapt_resp(im,:,1),prefTestOnly_noadapt_resp(im,:,2),'-o')
    hold on
    errorbar(testCons, prefTestOnly_contadapt_resp(im,:,1),prefTestOnly_contadapt_resp(im,:,2),'-o')
    errorbar(testCons, prefTestOnly_asynadapt_resp(im,:,1),prefTestOnly_asynadapt_resp(im,:,2),'-o')
    title(['Mask = ' num2str(maskCons(im))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
end
legend('No adapt', 'Contingent', 'Asynchronous');
suptitle(['Test only preferring cells- n = ' num2str(length([prefTestInd]))])
print(fullfile(fnout, 'adaptSummary_byMask_prefTestAll_250ms.pdf'),'-dpdf','-bestfit');


figure;
for im = 1:nMask
    subplot(2,3,im)
    errorbar(testCons, prefMaskOnly_noadapt_resp(:,im,1),prefMaskOnly_noadapt_resp(:,im,2),'-o')
    hold on
    errorbar(testCons, prefMaskOnly_contadapt_resp(:,im,1),prefMaskOnly_contadapt_resp(:,im,2),'-o')
    errorbar(testCons, prefMaskOnly_asynadapt_resp(:,im,1),prefMaskOnly_asynadapt_resp(:,im,2),'-o')
    title(['Test = ' num2str(maskCons(im))])
    ylabel('dF/F')
    xlabel('Mask Contrast')
    ylim([-0.1 0.5])
end
legend('No adapt', 'Contingent', 'Asynchronous');
suptitle(['Mask only preferring cells- n = ' num2str(length([prefMaskInd]))])
print(fullfile(fnout, 'adaptSummary_byMask_prefMaskAll_250ms.pdf'),'-dpdf','-bestfit');

figure;
for im = 1:nMask
    subplot(2,3,im)
    errorbar(testCons, prefPlaidOnly_noadapt_resp(im,:,1),prefPlaidOnly_noadapt_resp(im,:,2),'-o')
    hold on
    errorbar(testCons, prefPlaidOnly_contadapt_resp(im,:,1),prefPlaidOnly_contadapt_resp(im,:,2),'-o')
    errorbar(testCons, prefPlaidOnly_asynadapt_resp(im,:,1),prefPlaidOnly_asynadapt_resp(im,:,2),'-o')
    title(['Mask = ' num2str(maskCons(im))])
    ylabel('dF/F')
    xlabel('Test Contrast')
    ylim([-0.1 0.5])
end
legend('No adapt', 'Contingent', 'Asynchronous');
suptitle(['Plaid only preferring cells- n = ' num2str(length([prefPlaidInd]))])
print(fullfile(fnout, 'adaptSummary_byMask_prefPlaidOnly_250ms.pdf'),'-dpdf','-bestfit');

%%
prefTestOnly_noadapt_respCells = zeros(length(prefTestInd),nMask,nTest);
prefMaskOnly_noadapt_respCells = zeros(length(prefMaskInd),nMask,nTest);
prefPlaidOnly_noadapt_respCells = zeros(length(prefPlaidInd),nMask,nTest);
prefTestOrMask_noadapt_respCells = zeros(length(prefTestInd)+length(prefMaskInd),nMask,nTest);
prefTestOnly_contadapt_respCells = zeros(length(prefTestInd),nMask,nTest);
prefMaskOnly_contadapt_respCells = zeros(length(prefMaskInd),nMask,nTest);
prefPlaidOnly_contadapt_respCells = zeros(length(prefPlaidInd),nMask,nTest);
prefTestOrMask_contadapt_respCells = zeros(length(prefTestInd)+length(prefMaskInd),nMask,nTest);
prefTestOnly_asynadapt_respCells = zeros(length(prefTestInd),nMask,nTest);
prefMaskOnly_asynadapt_respCells = zeros(length(prefMaskInd),nMask,nTest);
prefPlaidOnly_asynadapt_respCells = zeros(length(prefPlaidInd),nMask,nTest);
prefTestOrMask_asynadapt_respCells = zeros(length(prefTestInd)+length(prefMaskInd),nMask,nTest);
for im = 1:nMask
    for it = 1:nTest
        prefTestOnly_noadapt_respCells(:,im,it) = noadapt_resp_cell_all{im,it}(prefTestInd,:)-noadapt_resp_cell_all{1,1}(prefTestInd,:);
        prefMaskOnly_noadapt_respCells(:,im,it) = noadapt_resp_cell_all{it,im}(prefMaskInd,:)-noadapt_resp_cell_all{1,1}(prefMaskInd,:);
        prefPlaidOnly_noadapt_respCells(:,im,it) = noadapt_resp_cell_all{im,it}(prefPlaidInd,:)-noadapt_resp_cell_all{1,1}(prefPlaidInd,:);
        prefTestOnly_contadapt_respCells(:,im,it) = contadapt_resp_cell_all{im,it}(prefTestInd,:)-contadapt_resp_cell_all{1,1}(prefTestInd,:);
        prefMaskOnly_contadapt_respCells(:,im,it) = contadapt_resp_cell_all{it,im}(prefMaskInd,:)-contadapt_resp_cell_all{1,1}(prefMaskInd,:);
        prefPlaidOnly_contadapt_respCells(:,im,it) = contadapt_resp_cell_all{im,it}(prefPlaidInd,:)-contadapt_resp_cell_all{1,1}(prefPlaidInd,:);
        prefTestOnly_asynadapt_respCells(:,im,it) = asynadapt_resp_cell_all{im,it}(prefTestInd,:)-asynadapt_resp_cell_all{1,1}(prefTestInd,:);
        prefMaskOnly_asynadapt_respCells(:,im,it) = asynadapt_resp_cell_all{it,im}(prefMaskInd,:)-asynadapt_resp_cell_all{1,1}(prefMaskInd,:);
        prefPlaidOnly_asynadapt_respCells(:,im,it) = asynadapt_resp_cell_all{im,it}(prefPlaidInd,:)-asynadapt_resp_cell_all{1,1}(prefPlaidInd,:);
        
        prefTestOrMask_noadapt_respCells(:,im,it) = [noadapt_resp_cell_all{im,it}(prefTestInd,:)-noadapt_resp_cell_all{1,1}(prefTestInd,:); noadapt_resp_cell_all{it,im}(prefMaskInd,:)-noadapt_resp_cell_all{1,1}(prefMaskInd,:)];
        prefTestOrMask_contadapt_respCells(:,im,it) = [contadapt_resp_cell_all{im,it}(prefTestInd,:)-contadapt_resp_cell_all{1,1}(prefTestInd,:); contadapt_resp_cell_all{it,im}(prefMaskInd,:)-contadapt_resp_cell_all{1,1}(prefMaskInd,:)];
        prefTestOrMask_asynadapt_respCells(:,im,it) = [asynadapt_resp_cell_all{im,it}(prefTestInd,:)-asynadapt_resp_cell_all{1,1}(prefTestInd,:); asynadapt_resp_cell_all{it,im}(prefMaskInd,:)-asynadapt_resp_cell_all{1,1}(prefMaskInd,:)];
    end
end

figure; 
subplot(1,2,1)
imagesc(squeeze(mean(prefTestOrMask_contadapt_respCells,1)))
title('Contingent')
colormap redblue
clim([-0.2 0.2])
axis square
xlabel('Test Contrast')
xticklabels(chop(testCons,2))
ylabel('Mask Contrast')
yticklabels(chop(maskCons,2))
colorbar
subplot(1,2,2)
imagesc(squeeze(mean(prefTestOrMask_asynadapt_respCells,1)))
title('Asynchronous')
colormap redblue
clim([-0.2 0.2])
axis square
xlabel('Test Contrast')
xticklabels(chop(testCons,2))
ylabel('Mask Contrast')
yticklabels(chop(maskCons,2))
colorbar
suptitle(['Test or Mask preferring cells- n = ' num2str(length([prefTestInd; prefMaskInd]))])
print(fullfile(fnout, 'adaptSummary_avgRespHeatmap_prefTestOrMaskAll_250ms.pdf'),'-dpdf','-bestfit');

% %subtract minimum (instead of no stim condition)
% for i = 1:size(prefTestOnly_noadapt_respCells,1)
%     prefTestOnly_noadapt_respCells(i,:,:) = prefTestOnly_noadapt_respCells(i,:,:)-min(min(prefTestOnly_noadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefMaskOnly_noadapt_respCells,1)
%     prefMaskOnly_noadapt_respCells(i,:,:) = prefMaskOnly_noadapt_respCells(i,:,:)-min(min(prefMaskOnly_noadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefPlaidOnly_noadapt_respCells,1)
%     prefPlaidOnly_noadapt_respCells(i,:,:) = prefPlaidOnly_noadapt_respCells(i,:,:)-min(min(prefPlaidOnly_noadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefTestOnly_contadapt_respCells,1)
%     prefTestOnly_contadapt_respCells(i,:,:) = prefTestOnly_contadapt_respCells(i,:,:)-min(min(prefTestOnly_contadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefMaskOnly_contadapt_respCells,1)
%     prefMaskOnly_contadapt_respCells(i,:,:) = prefMaskOnly_contadapt_respCells(i,:,:)-min(min(prefMaskOnly_contadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefPlaidOnly_contadapt_respCells,1)
%     prefPlaidOnly_contadapt_respCells(i,:,:) = prefPlaidOnly_contadapt_respCells(i,:,:)-min(min(prefPlaidOnly_contadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefTestOnly_asynadapt_respCells,1)
%     prefTestOnly_asynadapt_respCells(i,:,:) = prefTestOnly_asynadapt_respCells(i,:,:)-min(min(prefTestOnly_asynadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefMaskOnly_asynadapt_respCells,1)
%     prefMaskOnly_asynadapt_respCells(i,:,:) = prefMaskOnly_asynadapt_respCells(i,:,:)-min(min(prefMaskOnly_asynadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefPlaidOnly_asynadapt_respCells,1)
%     prefPlaidOnly_asynadapt_respCells(i,:,:) = prefPlaidOnly_asynadapt_respCells(i,:,:)-min(min(prefPlaidOnly_asynadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefTestOrMask_noadapt_respCells,1)
%     prefTestOrMask_noadapt_respCells(i,:,:) = prefTestOrMask_noadapt_respCells(i,:,:)-min(min(prefTestOrMask_noadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefTestOrMask_contadapt_respCells,1)
%     prefTestOrMask_contadapt_respCells(i,:,:) = prefTestOrMask_contadapt_respCells(i,:,:)-min(min(prefTestOrMask_contadapt_respCells(i,:,:)));
% end
% for i = 1:size(prefTestOrMask_asynadapt_respCells,1)
%     prefTestOrMask_asynadapt_respCells(i,:,:) = prefTestOrMask_asynadapt_respCells(i,:,:)-min(min(prefTestOrMask_asynadapt_respCells(i,:,:)));
% end

% %subtract mask-alone response
% for im = 1:nMask
%     prefTestOnly_noadapt_respCells(:,im,:) = prefTestOnly_noadapt_respCells(:,im,:)-prefTestOnly_noadapt_respCells(:,im,1);
%     prefMaskOnly_noadapt_respCells(:,im,:) = prefMaskOnly_noadapt_respCells(:,im,:)-prefMaskOnly_noadapt_respCells(:,im,1);
%     prefPlaidOnly_noadapt_respCells(:,im,:) =  prefPlaidOnly_noadapt_respCells(:,im,:)-prefPlaidOnly_noadapt_respCells(:,im,1);
%     prefTestOnly_contadapt_respCells(:,im,:) = prefTestOnly_contadapt_respCells(:,im,:)-prefTestOnly_contadapt_respCells(:,im,1);
%     prefMaskOnly_contadapt_respCells(:,im,:) =  prefMaskOnly_contadapt_respCells(:,im,:)- prefMaskOnly_contadapt_respCells(:,im,1);
%     prefPlaidOnly_contadapt_respCells(:,im,:) = prefPlaidOnly_contadapt_respCells(:,im,:)-prefPlaidOnly_contadapt_respCells(:,im,1);
%     prefTestOnly_asynadapt_respCells(:,im,:) = prefTestOnly_asynadapt_respCells(:,im,:)-prefTestOnly_asynadapt_respCells(:,im,1);
%     prefMaskOnly_asynadapt_respCells(:,im,:) =  prefMaskOnly_asynadapt_respCells(:,im,:)- prefMaskOnly_asynadapt_respCells(:,im,1);
%     prefPlaidOnly_asynadapt_respCells(:,im,:) = prefPlaidOnly_asynadapt_respCells(:,im,:)-prefPlaidOnly_asynadapt_respCells(:,im,1);
%     prefTestOrMask_noadapt_respCells(:,im,:) = prefTestOrMask_noadapt_respCells(:,im,:)-prefTestOrMask_noadapt_respCells(:,im,1);
%     prefTestOrMask_contadapt_respCells(:,im,:) = prefTestOrMask_contadapt_respCells(:,im,:)-prefTestOrMask_contadapt_respCells(:,im,1);
%     prefTestOrMask_asynadapt_respCells(:,im,:) = prefTestOrMask_asynadapt_respCells(:,im,:)-prefTestOrMask_asynadapt_respCells(:,im,1);
% end

%set negative values to 0
prefTestOnly_noadapt_respCells(find(prefTestOnly_noadapt_respCells<0)) = 0;
prefMaskOnly_noadapt_respCells(find(prefMaskOnly_noadapt_respCells<0)) = 0;
prefPlaidOnly_noadapt_respCells(find(prefPlaidOnly_noadapt_respCells<0)) = 0;
prefTestOnly_contadapt_respCells(find(prefTestOnly_contadapt_respCells<0)) = 0;
prefMaskOnly_contadapt_respCells(find(prefMaskOnly_contadapt_respCells<0)) = 0;
prefPlaidOnly_contadapt_respCells(find(prefPlaidOnly_contadapt_respCells<0)) = 0;
prefTestOnly_asynadapt_respCells(find(prefTestOnly_asynadapt_respCells<0)) = 0;
prefMaskOnly_asynadapt_respCells(find(prefMaskOnly_asynadapt_respCells<0)) = 0;
prefPlaidOnly_asynadapt_respCells(find(prefPlaidOnly_asynadapt_respCells<0)) = 0;
prefTestOrMask_noadapt_respCells(find(prefTestOrMask_noadapt_respCells<0)) = 0;
prefTestOrMask_contadapt_respCells(find(prefTestOrMask_contadapt_respCells<0)) = 0;
prefTestOrMask_asynadapt_respCells(find(prefTestOrMask_asynadapt_respCells<0)) = 0;



        
prefTestOnly_noadapt_MI = zeros(length(prefTestInd),nTest,nMask-1);
prefMaskOnly_noadapt_MI = zeros(length(prefMaskInd),nTest,nMask-1);
prefPlaidOnly_noadapt_MI = zeros(length(prefPlaidInd),nTest,nMask-1);
prefTestOrMask_noadapt_MI = zeros(length(prefTestInd)+length(prefMaskInd),nTest,nMask-1);
prefTestOnly_asynadapt_MI = zeros(length(prefTestInd),nTest,nMask-1);
prefMaskOnly_asynadapt_MI = zeros(length(prefMaskInd),nTest,nMask-1);
prefPlaidOnly_asynadapt_MI = zeros(length(prefPlaidInd),nTest,nMask-1);
prefTestOrMask_asynadapt_MI = zeros(length(prefTestInd)+length(prefMaskInd),nTest,nMask-1);
prefTestOnly_contadapt_MI = zeros(length(prefTestInd),nTest,nMask-1);
prefMaskOnly_contadapt_MI = zeros(length(prefMaskInd),nTest,nMask-1);
prefPlaidOnly_contadapt_MI = zeros(length(prefPlaidInd),nTest,nMask-1);
prefTestOrMask_contadapt_MI = zeros(length(prefTestInd)+length(prefMaskInd),nTest,nMask-1);
    
for im = 2:nMask
    prefTestOnly_noadapt_MI(:,:,im-1) = squeeze(prefTestOnly_noadapt_respCells(:,1,:)-prefTestOnly_noadapt_respCells(:,im,:))./squeeze(prefTestOnly_noadapt_respCells(:,1,:)+prefTestOnly_noadapt_respCells(:,im,:));
    prefMaskOnly_noadapt_MI(:,:,im-1) = squeeze(prefMaskOnly_noadapt_respCells(:,1,:)-prefMaskOnly_noadapt_respCells(:,im,:))./squeeze(prefMaskOnly_noadapt_respCells(:,1,:)+prefMaskOnly_noadapt_respCells(:,im,:));
    prefPlaidOnly_noadapt_MI(:,:,im-1) = squeeze(prefPlaidOnly_noadapt_respCells(:,1,:)-prefPlaidOnly_noadapt_respCells(:,im,:))./squeeze(prefPlaidOnly_noadapt_respCells(:,1,:)+prefPlaidOnly_noadapt_respCells(:,im,:));
    prefTestOrMask_noadapt_MI(:,:,im-1) = squeeze(prefTestOrMask_noadapt_respCells(:,1,:)-prefTestOrMask_noadapt_respCells(:,im,:))./squeeze(prefTestOrMask_noadapt_respCells(:,1,:)+prefTestOrMask_noadapt_respCells(:,im,:));
    
    prefTestOnly_asynadapt_MI(:,:,im-1) = squeeze(prefTestOnly_asynadapt_respCells(:,1,:)-prefTestOnly_asynadapt_respCells(:,im,:))./squeeze(prefTestOnly_asynadapt_respCells(:,1,:)+prefTestOnly_asynadapt_respCells(:,im,:));
    prefMaskOnly_asynadapt_MI(:,:,im-1) = squeeze(prefMaskOnly_asynadapt_respCells(:,1,:)-prefMaskOnly_asynadapt_respCells(:,im,:))./squeeze(prefMaskOnly_asynadapt_respCells(:,1,:)+prefMaskOnly_asynadapt_respCells(:,im,:));
    prefPlaidOnly_asynadapt_MI(:,:,im-1) = squeeze(prefPlaidOnly_asynadapt_respCells(:,1,:)-prefPlaidOnly_asynadapt_respCells(:,im,:))./squeeze(prefPlaidOnly_asynadapt_respCells(:,1,:)+prefPlaidOnly_asynadapt_respCells(:,im,:));
    prefTestOrMask_asynadapt_MI(:,:,im-1) = squeeze(prefTestOrMask_asynadapt_respCells(:,1,:)-prefTestOrMask_asynadapt_respCells(:,im,:))./squeeze(prefTestOrMask_asynadapt_respCells(:,1,:)+prefTestOrMask_asynadapt_respCells(:,im,:));
    
    prefTestOnly_contadapt_MI(:,:,im-1) = squeeze(prefTestOnly_contadapt_respCells(:,1,:)-prefTestOnly_contadapt_respCells(:,im,:))./squeeze(prefTestOnly_contadapt_respCells(:,1,:)+prefTestOnly_contadapt_respCells(:,im,:));
    prefMaskOnly_contadapt_MI(:,:,im-1) = squeeze(prefMaskOnly_contadapt_respCells(:,1,:)-prefMaskOnly_contadapt_respCells(:,im,:))./squeeze(prefMaskOnly_contadapt_respCells(:,1,:)+prefMaskOnly_contadapt_respCells(:,im,:));
    prefPlaidOnly_contadapt_MI(:,:,im-1) = squeeze(prefPlaidOnly_contadapt_respCells(:,1,:)-prefPlaidOnly_contadapt_respCells(:,im,:))./squeeze(prefPlaidOnly_contadapt_respCells(:,1,:)+prefPlaidOnly_contadapt_respCells(:,im,:));
    prefTestOrMask_contadapt_MI(:,:,im-1) = squeeze(prefTestOrMask_contadapt_respCells(:,1,:)-prefTestOrMask_contadapt_respCells(:,im,:))./squeeze(prefTestOrMask_contadapt_respCells(:,1,:)+prefTestOrMask_contadapt_respCells(:,im,:));
end


%Test pref
figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefTestOnly_noadapt_MI(:,it,:)),1);
    n_asynadapt = sum(~isnan(prefTestOnly_asynadapt_MI(:,it,:)),1);
    n_contadapt = sum(~isnan(prefTestOnly_asynadapt_MI(:,it,:)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOnly_noadapt_MI(:,it,:),1)), squeeze(nanstd(prefTestOnly_noadapt_MI(:,it,:),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOnly_asynadapt_MI(:,it,:),1)), squeeze(nanstd(prefTestOnly_asynadapt_MI(:,it,:),1)./sqrt(n_asynadapt)), '-o')
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOnly_contadapt_MI(:,it,:),1)), squeeze(nanstd(prefTestOnly_contadapt_MI(:,it,:),1)./sqrt(n_contadapt)), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('MI')
    ylim([-1 1])
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefTestOnly_noadapt_MI,2)),1);
n_asynadapt = sum(~isnan(nanmean(prefTestOnly_asynadapt_MI,2)),1);
n_contadapt = sum(~isnan(nanmean(prefTestOnly_contadapt_MI,2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOnly_noadapt_MI,2),1)), squeeze(nanstd(nanmean(prefTestOnly_noadapt_MI,2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOnly_asynadapt_MI,2),1)), squeeze(nanstd(nanmean(prefTestOnly_asynadapt_MI,2),1)./sqrt(n_asynadapt)), '-o')
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOnly_contadapt_MI,2),1)), squeeze(nanstd(nanmean(prefTestOnly_contadapt_MI,2),1)./sqrt(n_contadapt)), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('MI')
ylim([-1 1])
subplot(3,2,6)
n_asynadapt = squeeze(sum(~isnan(nanmean(prefTestOnly_asynadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2)),1));
n_contadapt = squeeze(sum(~isnan(nanmean(prefTestOnly_contadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefTestOnly_asynadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefTestOnly_asynadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2),[],1))./sqrt(n_asynadapt)),'-o');
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefTestOnly_contadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefTestOnly_contadapt_MI,2)-nanmean(prefTestOnly_noadapt_MI,2),[],1))./sqrt(n_contadapt)),'-o');
title('All Test')
xlabel('Mask contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
legend({'No adapt', 'Asynchronous', 'Contingent'})
suptitle(['Test Pref- n = ' num2str(min(n_asynadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_MI_prefTestAll_250ms.pdf'),'-dpdf','-bestfit');
%Mask pref
figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefMaskOnly_noadapt_MI(:,it,:)),1);
    n_asynadapt = sum(~isnan(prefMaskOnly_asynadapt_MI(:,it,:)),1);
    n_contadapt = sum(~isnan(prefMaskOnly_asynadapt_MI(:,it,:)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefMaskOnly_noadapt_MI(:,it,:),1)), squeeze(nanstd(prefMaskOnly_noadapt_MI(:,it,:),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefMaskOnly_asynadapt_MI(:,it,:),1)), squeeze(nanstd(prefMaskOnly_asynadapt_MI(:,it,:),1)./sqrt(n_asynadapt)), '-o')
    errorbar(maskCons(2:5), squeeze(nanmean(prefMaskOnly_contadapt_MI(:,it,:),1)), squeeze(nanstd(prefMaskOnly_contadapt_MI(:,it,:),1)./sqrt(n_contadapt)), '-o')
    title(['Mask = ' num2str(testCons(it))])
    xlabel('Test contrast')
    ylabel('MI')
    ylim([-1 1])
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefMaskOnly_noadapt_MI,2)),1);
n_asynadapt = sum(~isnan(nanmean(prefMaskOnly_asynadapt_MI,2)),1);
n_contadapt = sum(~isnan(nanmean(prefMaskOnly_contadapt_MI,2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefMaskOnly_noadapt_MI,2),1)), squeeze(nanstd(nanmean(prefMaskOnly_noadapt_MI,2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefMaskOnly_asynadapt_MI,2),1)), squeeze(nanstd(nanmean(prefMaskOnly_asynadapt_MI,2),1)./sqrt(n_asynadapt)), '-o')
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefMaskOnly_contadapt_MI,2),1)), squeeze(nanstd(nanmean(prefMaskOnly_contadapt_MI,2),1)./sqrt(n_contadapt)), '-o')
title('All Mask')
xlabel('Test contrast')
ylabel('MI')
ylim([-1 1])
subplot(3,2,6)
n_asynadapt = squeeze(sum(~isnan(nanmean(prefMaskOnly_asynadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2)),1));
n_contadapt = squeeze(sum(~isnan(nanmean(prefMaskOnly_contadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefMaskOnly_asynadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefMaskOnly_asynadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2),[],1))./sqrt(n_asynadapt)),'-o');
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefMaskOnly_contadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefMaskOnly_contadapt_MI,2)-nanmean(prefMaskOnly_noadapt_MI,2),[],1))./sqrt(n_contadapt)),'-o');
title('All Mask')
xlabel('Test contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
legend({'No adapt', 'Asynchronous', 'Contingent'})
suptitle(['Mask Pref- n = ' num2str(min(n_asynadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_MI_prefMaskAll_250ms.pdf'),'-dpdf','-bestfit');

%Test or Mask pref
figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefTestOrMask_noadapt_MI(:,it,:)),1);
    n_asynadapt = sum(~isnan(prefTestOrMask_asynadapt_MI(:,it,:)),1);
    n_contadapt = sum(~isnan(prefTestOrMask_asynadapt_MI(:,it,:)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOrMask_noadapt_MI(:,it,:),1)), squeeze(nanstd(prefTestOrMask_noadapt_MI(:,it,:),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOrMask_asynadapt_MI(:,it,:),1)), squeeze(nanstd(prefTestOrMask_asynadapt_MI(:,it,:),1)./sqrt(n_asynadapt)), '-o')
    errorbar(maskCons(2:5), squeeze(nanmean(prefTestOrMask_contadapt_MI(:,it,:),1)), squeeze(nanstd(prefTestOrMask_contadapt_MI(:,it,:),1)./sqrt(n_contadapt)), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('MI')
    ylim([-1 1])
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefTestOrMask_noadapt_MI,2)),1);
n_asynadapt = sum(~isnan(nanmean(prefTestOrMask_asynadapt_MI,2)),1);
n_contadapt = sum(~isnan(nanmean(prefTestOrMask_contadapt_MI,2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOrMask_noadapt_MI,2),1)), squeeze(nanstd(nanmean(prefTestOrMask_noadapt_MI,2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOrMask_asynadapt_MI,2),1)), squeeze(nanstd(nanmean(prefTestOrMask_asynadapt_MI,2),1)./sqrt(n_asynadapt)), '-o')
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefTestOrMask_contadapt_MI,2),1)), squeeze(nanstd(nanmean(prefTestOrMask_contadapt_MI,2),1)./sqrt(n_contadapt)), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('MI')
ylim([-1 1])
subplot(3,2,6)
n_asynadapt = squeeze(sum(~isnan(nanmean(prefTestOrMask_asynadapt_MI,2)-nanmean(prefTestOrMask_noadapt_MI,2)),1));
n_contadapt = squeeze(sum(~isnan(nanmean(prefTestOrMask_contadapt_MI,2)-nanmean(prefTestOrMask_noadapt_MI,2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefTestOrMask_asynadapt_MI,2)-nanmean(prefTestOrMask_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefTestOrMask_asynadapt_MI,2)-nanmean(prefTestOrMask_noadapt_MI,2),[],1))./sqrt(n_asynadapt)),'-o');
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefTestOrMask_contadapt_MI,2)-nanmean(prefTestOrMask_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefTestOrMask_contadapt_MI,2)-nanmean(prefTestOrMask_noadapt_MI,2),[],1))./sqrt(n_contadapt)),'-o');
title('All Test')
xlabel('Mask contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
legend({'No adapt', 'Asynchronous', 'Contingent'})
suptitle(['Test or Mask Pref- n = ' num2str(min(n_asynadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_MI_prefTestOrMaskAll_250ms.pdf'),'-dpdf','-bestfit');

figure;
for it = 2:nMask
    subplot(3,2,it-1)
    n_noadapt = sum(~isnan(prefPlaidOnly_noadapt_MI(:,it,:)),1);
    n_asynadapt = sum(~isnan(prefPlaidOnly_asynadapt_MI(:,it,:)),1);
    n_contadapt = sum(~isnan(prefPlaidOnly_asynadapt_MI(:,it,:)),1);
    errorbar(maskCons(2:5), squeeze(nanmean(prefPlaidOnly_noadapt_MI(:,it,:),1)), squeeze(nanstd(prefPlaidOnly_noadapt_MI(:,it,:),1)./sqrt(n_noadapt)), '-o')
    hold on
    errorbar(maskCons(2:5), squeeze(nanmean(prefPlaidOnly_asynadapt_MI(:,it,:),1)), squeeze(nanstd(prefPlaidOnly_asynadapt_MI(:,it,:),1)./sqrt(n_asynadapt)), '-o')
    errorbar(maskCons(2:5), squeeze(nanmean(prefPlaidOnly_contadapt_MI(:,it,:),1)), squeeze(nanstd(prefPlaidOnly_contadapt_MI(:,it,:),1)./sqrt(n_contadapt)), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('MI')
    ylim([-1 1])
end
subplot(3,2,5)
n_noadapt = sum(~isnan(nanmean(prefPlaidOnly_noadapt_MI,2)),1);
n_asynadapt = sum(~isnan(nanmean(prefPlaidOnly_asynadapt_MI,2)),1);
n_contadapt = sum(~isnan(nanmean(prefPlaidOnly_contadapt_MI,2)),1);
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefPlaidOnly_noadapt_MI,2),1)), squeeze(nanstd(nanmean(prefPlaidOnly_noadapt_MI,2),1)./sqrt(n_noadapt)), '-o')
hold on
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefPlaidOnly_asynadapt_MI,2),1)), squeeze(nanstd(nanmean(prefPlaidOnly_asynadapt_MI,2),1)./sqrt(n_asynadapt)), '-o')
errorbar(maskCons(2:5), squeeze(nanmean(nanmean(prefPlaidOnly_contadapt_MI,2),1)), squeeze(nanstd(nanmean(prefPlaidOnly_contadapt_MI,2),1)./sqrt(n_contadapt)), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('MI')
ylim([-1 1])
subplot(3,2,6)
n_asynadapt = squeeze(sum(~isnan(nanmean(prefPlaidOnly_asynadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2)),1));
n_contadapt = squeeze(sum(~isnan(nanmean(prefPlaidOnly_contadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2)),1));
plot(nan,nan)
hold on
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefPlaidOnly_asynadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefPlaidOnly_asynadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2),[],1))./sqrt(n_asynadapt)),'-o');
errorbar(maskCons(2:5),squeeze(nanmean(nanmean(prefPlaidOnly_contadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2),1)), (squeeze(nanstd(nanmean(prefPlaidOnly_contadapt_MI,2)-nanmean(prefPlaidOnly_noadapt_MI,2),[],1))./sqrt(n_contadapt)),'-o');
title('All Test')
xlabel('Mask contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
legend({'No adapt', 'Asynchronous', 'Contingent'})
suptitle(['Plaid Pref- n = ' num2str(min(n_asynadapt,[],1))])
print(fullfile(fnout, 'adaptSummary_MI_prefPlaidOnly_250ms.pdf'),'-dpdf','-bestfit');

%%
prefTestOrMask_noadapt_respdiag = zeros(nTest,2);
prefTestOrMask_noadapt_respTOnly = zeros(nTest,2);
prefPlaidOnly_noadapt_respdiag = zeros(nTest,2);
prefPlaidOnly_noadapt_respTOnly = zeros(nTest,2);
prefTestOrMask_contadapt_respdiag = zeros(nTest,2);
prefTestOrMask_contadapt_respTOnly = zeros(nTest,2);
prefPlaidOnly_contadapt_respdiag = zeros(nTest,2);
prefPlaidOnly_contadapt_respTOnly = zeros(nTest,2);
prefTestOrMask_asynadapt_respdiag = zeros(nTest,2);
prefTestOrMask_asynadapt_respTOnly = zeros(nTest,2);
prefPlaidOnly_asynadapt_respdiag = zeros(nTest,2);
prefPlaidOnly_asynadapt_respTOnly = zeros(nTest,2);
for it = 1:nTest
    prefTestOrMask_noadapt_respdiag(it,:) = prefTestOrMask_noadapt_resp(it,it,:);
    prefTestOrMask_noadapt_respTOnly(it,:) = prefTestOrMask_noadapt_resp(1,it,:);
    prefPlaidOnly_noadapt_respdiag(it,:) = prefPlaidOnly_noadapt_resp(it,it,:);
    prefPlaidOnly_noadapt_respTOnly(it,:) = prefPlaidOnly_noadapt_resp(1,it,:);
    prefTestOrMask_contadapt_respdiag(it,:) = prefTestOrMask_contadapt_resp(it,it,:);
    prefTestOrMask_contadapt_respTOnly(it,:) = prefTestOrMask_contadapt_resp(1,it,:);
    prefPlaidOnly_contadapt_respdiag(it,:) = prefPlaidOnly_contadapt_resp(it,it,:);
    prefPlaidOnly_contadapt_respTOnly(it,:) = prefPlaidOnly_contadapt_resp(1,it,:);
    prefTestOrMask_asynadapt_respdiag(it,:) = prefTestOrMask_asynadapt_resp(it,it,:);
    prefTestOrMask_asynadapt_respTOnly(it,:) = prefTestOrMask_asynadapt_resp(1,it,:);
    prefPlaidOnly_asynadapt_respdiag(it,:) = prefPlaidOnly_asynadapt_resp(it,it,:);
    prefPlaidOnly_asynadapt_respTOnly(it,:) = prefPlaidOnly_asynadapt_resp(1,it,:);
end
figure;
subplot(1,3,1)
errorbar(maskCons,prefTestOrMask_noadapt_respdiag(:,1),prefTestOrMask_noadapt_respdiag(:,2),'-o')
hold on
errorbar(maskCons,prefTestOrMask_noadapt_respTOnly(:,1),prefTestOrMask_noadapt_respTOnly(:,2),'-o')
errorbar(maskCons,prefPlaidOnly_noadapt_respdiag(:,1),prefPlaidOnly_noadapt_respdiag(:,2),'-o')
errorbar(maskCons,prefPlaidOnly_noadapt_respTOnly(:,1),prefPlaidOnly_noadapt_respTOnly(:,2),'-o')
title('No adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,2)
errorbar(maskCons,prefTestOrMask_contadapt_respdiag(:,1),prefTestOrMask_contadapt_respdiag(:,2),'-o')
hold on
errorbar(maskCons,prefTestOrMask_contadapt_respTOnly(:,1),prefTestOrMask_contadapt_respTOnly(:,2),'-o')
errorbar(maskCons,prefPlaidOnly_contadapt_respdiag(:,1),prefPlaidOnly_contadapt_respdiag(:,2),'-o')
errorbar(maskCons,prefPlaidOnly_contadapt_respTOnly(:,1),prefPlaidOnly_contadapt_respTOnly(:,2),'-o')
title('Contingent adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
subplot(1,3,3)
errorbar(maskCons,prefTestOrMask_asynadapt_respdiag(:,1),prefTestOrMask_asynadapt_respdiag(:,2),'-o')
hold on
errorbar(maskCons,prefTestOrMask_asynadapt_respTOnly(:,1),prefTestOrMask_asynadapt_respTOnly(:,2),'-o')
errorbar(maskCons,prefPlaidOnly_asynadapt_respdiag(:,1),prefPlaidOnly_asynadapt_respdiag(:,2),'-o')
errorbar(maskCons,prefPlaidOnly_asynadapt_respTOnly(:,1),prefPlaidOnly_asynadapt_respTOnly(:,2),'-o')
title('Asynchronous adapt')
ylabel('dF/F')
xlabel('Test Contrast')
ylim([-0.1 0.5])
legend({'Test/Mask- diag','Test/Mask- TestOnly','Plaid- diag', 'Plaid- Test Only'})
suptitle(['Diagonal vs TestOnly stim'])
print(fullfile(fnout,'Summary_DaigvTestAll_250ms.pdf'),'-dpdf','-bestfit');

%% ori tuning and supp index
resp_single = unique([resptest_ind_all; respmask_ind_all]);
resp_any = unique([resptest_ind_all; respmask_ind_all; respplaid_ind_all]);
resp_plaid_only = setdiff(respplaid_ind_all,resp_single);
noadapt_rect_resp = noadapt_resp_cell_all;
for im = 1:nMask
    for it = 1:nTest
        noadapt_rect_resp{im,it}(find(noadapt_rect_resp{im,it}<0)) = 0;
    end
end
resp_all = [noadapt_rect_resp{1,5} noadapt_rect_resp{5,1} noadapt_rect_resp{5,5}];
ind0 = find(noadapt_rect_resp{1,5}+noadapt_rect_resp{5,1}==0);
resp_single = setdiff(resp_single,ind0);
sum_all = noadapt_rect_resp{1,5}+noadapt_rect_resp{5,1};
osi_all = abs(noadapt_rect_resp{1,5}-noadapt_rect_resp{5,1})./(noadapt_rect_resp{1,5}+noadapt_rect_resp{5,1});
osi_all(ind0) = 0;
supp_all = noadapt_rect_resp{5,5}./(noadapt_rect_resp{1,5}+noadapt_rect_resp{5,1});
resp_thresh_ind = unique([find(noadapt_rect_resp{1,5}>0.05); find(noadapt_rect_resp{5,5}>0.05); find(noadapt_rect_resp{5,1}>0.05)]);

ind_facil = find(supp_all(resp_single)>1);
ind_supp = find(supp_all(resp_single)<=1);
figure; 
subplot(2,2,1)
scatter(osi_all(resp_single(ind_supp)), supp_all(resp_single(ind_supp)));
xlabel('OSI')
ylabel('SI')
title(['SI <= 1 - n = ' num2str(length(ind_supp))])
subplot(2,2,2)
scatter(osi_all(resp_single(ind_facil)), supp_all(resp_single(ind_facil)));
set(gca, 'yscale','log')
xlabel('OSI')
ylabel('SI')
title(['SI > 1 - n = ' num2str(length(ind_facil))])
[n bin] = histc(osi_all,[0:0.15:1.2]);

for i = 1:length(n)
    ind = intersect(resp_single(ind_supp),find(bin==i));
    subplot(2,2,3)
    errorbarxy(mean(osi_all(ind),1), mean(supp_all(ind),1), std(osi_all(ind),[],1)./sqrt(length(ind)), std(supp_all(ind),[],1)./sqrt(length(ind)));
    hold on
    ind = intersect(resp_single(ind_facil),find(bin==i));    
    subplot(2,2,4)
    errorbarxy(mean(osi_all(ind),1), geomean(supp_all(ind),1), std(osi_all(ind),[],1)./sqrt(length(ind)), std(supp_all(ind),[],1)./sqrt(length(ind)));
    hold on
end
subplot(2,2,3)
xlabel('OSI')
ylabel('SI (mean)')
ylim([0 1])
xlim([0 1])
subplot(2,2,4)
xlabel('OSI')
ylabel('SI (geomean)')
xlim([0 1])
suptitle(['Single grating responsive cells- n = ' num2str(length(resp_single))]) 
print(fullfile(fnout,'Summary_NoAdapt_OSIvSI.pdf'),'-dpdf','-bestfit');

ind_facil_osi = intersect(find(supp_all(resp_single)>1), find(osi_all>0.7));
ind_supp_osi = intersect(find(supp_all(resp_single)<=1), find(osi_all>0.7));
rat = length(ind_facil_osi)./(length(ind_facil_osi)+length(ind_supp_osi));
