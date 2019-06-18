%contingent adaptation experiments
mouse_mat = strvcat('i1303','i1103','i1304');
date_mat = strvcat('190501','190508','190508');
run_mat = strvcat('002','002','002');

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
print(fullfile(fnout, 'adaptSummary_testPrefOnly.pdf'),'-dpdf','-bestfit');

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
print(fullfile(fnout, 'adaptSummary_maskPrefOnly.pdf'),'-dpdf','-bestfit');

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
print(fullfile(fnout, 'adaptSummary_plaidPrefOnly.pdf'),'-dpdf','-bestfit');

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
print(fullfile(fnout, 'adaptSummary_prefTestOrMask.pdf'),'-dpdf','-bestfit');

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
print(fullfile(fnout, 'adaptSummary_byMask_prefTestOrMask.pdf'),'-dpdf','-bestfit');

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
print(fullfile(fnout, 'adaptSummary_byMask_prefTestOnly.pdf'),'-dpdf','-bestfit');


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
print(fullfile(fnout, 'adaptSummary_byMask_prefMaskOnly.pdf'),'-dpdf','-bestfit');

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
print(fullfile(fnout, 'adaptSummary_byMask_prefPlaidOnly.pdf'),'-dpdf','-bestfit');

%%
T_noadapt_TestOrMask = prefTestOrMask_noadapt_resp(:,:,1);
T_noadapt_TestOrMask(find(T_noadapt_TestOrMask<0)) = 0;
MI_noadapt_TestOrMask = (T_noadapt_TestOrMask(1,:)-T_noadapt_TestOrMask(2:5,:))./(T_noadapt_TestOrMask(1,:)+T_noadapt_TestOrMask(2:5,:));

T_contadapt_TestOrMask = prefTestOrMask_contadapt_resp(:,:,1)-prefTestOrMask_contadapt_resp(5,1,1);
T_contadapt_TestOrMask(find(T_contadapt_TestOrMask<0)) = 0;
MI_contadapt_TestOrMask = (T_contadapt_TestOrMask(1,:)-T_contadapt_TestOrMask(2:5,:))./(T_contadapt_TestOrMask(1,:)+T_contadapt_TestOrMask(2:5,:));

T_asynadapt_TestOrMask = prefTestOrMask_asynadapt_resp(:,:,1)-prefTestOrMask_asynadapt_resp(5,1,1);
T_asynadapt_TestOrMask(find(T_asynadapt_TestOrMask<0)) = 0;
MI_asynadapt_TestOrMask = (T_asynadapt_TestOrMask(1,:)-T_asynadapt_TestOrMask(2:5,:))./(T_asynadapt_TestOrMask(1,:)+T_asynadapt_TestOrMask(2:5,:));

T_noadapt_Test = prefTestOnly_noadapt_resp(:,:,1);
T_noadapt_Test(find(T_noadapt_Test<0)) = 0;
MI_noadapt_Test = (T_noadapt_Test(1,:)-T_noadapt_Test(2:5,:))./(T_noadapt_Test(1,:)+T_noadapt_Test(2:5,:));

T_contadapt_Test = prefTestOnly_contadapt_resp(:,:,1)-prefTestOnly_contadapt_resp(5,1,1);
T_contadapt_Test(find(T_contadapt_Test<0)) = 0;
MI_contadapt_Test = (T_contadapt_Test(1,:)-T_contadapt_Test(2:5,:))./(T_contadapt_Test(1,:)+T_contadapt_Test(2:5,:));

T_asynadapt_Test = prefTestOnly_asynadapt_resp(:,:,1)-prefTestOnly_asynadapt_resp(5,1,1);
T_asynadapt_Test(find(T_asynadapt_Test<0)) = 0;
MI_asynadapt_Test = (T_asynadapt_Test(1,:)-T_asynadapt_Test(2:5,:))./(T_asynadapt_Test(1,:)+T_asynadapt_Test(2:5,:));

T_noadapt_Mask = prefMaskOnly_noadapt_resp(:,:,1);
T_noadapt_Mask(find(T_noadapt_Mask<0)) = 0;
MI_noadapt_Mask = (T_noadapt_Mask(:,1)-T_noadapt_Mask(:,2:5))./(T_noadapt_Mask(:,1)+T_noadapt_Mask(:,2:5));

T_contadapt_Mask = prefMaskOnly_contadapt_resp(:,:,1)-prefMaskOnly_contadapt_resp(1,5,1);
T_contadapt_Mask(find(T_contadapt_Mask<0)) = 0;
MI_contadapt_Mask = (T_contadapt_Mask(:,1)-T_contadapt_Mask(:,2:5))./(T_contadapt_Mask(:,1)+T_contadapt_Mask(:,2:5));

T_asynadapt_Mask = prefMaskOnly_asynadapt_resp(:,:,1)-prefMaskOnly_asynadapt_resp(1,5,1);
T_asynadapt_Mask(find(T_asynadapt_Mask<0)) = 0;
MI_asynadapt_Mask = (T_asynadapt_Mask(:,1)-T_asynadapt_Mask(:,2:5))./(T_asynadapt_Mask(:,1)+T_asynadapt_Mask(:,2:5));


% MI_cont_diff = nanmean(MI_contadapt_Mask-MI_noadapt_Mask,2);
% MI_asyn_diff = nanmean(MI_asynadapt_Mask-MI_noadapt_Mask,2);

figure; 
for it = 2:nMask
    subplot(3,2,it-1)
    plot(maskCons(2:5), MI_noadapt_TestOrMask(:,it), '-o')
    hold on
    plot(maskCons(2:5), MI_contadapt_TestOrMask(:,it), '-o')
    plot(maskCons(2:5), MI_asynadapt_TestOrMask(:,it), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('MI')
    ylim([-1 1])
end
subplot(3,2,5)
plot(maskCons(2:5),nanmean(MI_noadapt_TestOrMask(:,2:5),2),'-o')
hold on
plot(maskCons(2:5), nanmean(MI_contadapt_TestOrMask(:,2:5),2), '-o')
plot(maskCons(2:5), nanmean(MI_asynadapt_TestOrMask(:,2:5),2), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('MI')
ylim([-1 1])
subplot(3,2,6)
plot(nan,nan)
hold on
errorbar(maskCons(2:5),nanmean(MI_contadapt_TestOrMask(:,2:5)-MI_noadapt_TestOrMask(:,2:5),2),'-o')
errorbar(maskCons(2:5), nanmean(MI_asynadapt_TestOrMask(:,2:5)-MI_noadapt_TestOrMask(:,2:5),2), '-o')
title('All Test')
xlabel('Mask contrast')
ylabel('MI (post-pre)')
ylim([-1 1])
suptitle('Test or Mask Pref')

figure; 
for it = 2:nMask
    subplot(3,2,it-1)
    plot(maskCons(2:5), MI_noadapt_Test(:,it), '-o')
    hold on
    plot(maskCons(2:5), MI_contadapt_Test(:,it), '-o')
    plot(maskCons(2:5), MI_asynadapt_Test(:,it), '-o')
    title(['Test = ' num2str(testCons(it))])
    xlabel('Mask contrast')
    ylabel('MI')
end
suptitle('Test Pref')

figure; 
for it = 2:nMask
    subplot(3,2,it-1)
    plot(maskCons(2:5)', MI_noadapt_Mask(it,:), '-o')
    hold on
    plot(maskCons(2:5), MI_contadapt_Mask(it,:), '-o')
    plot(maskCons(2:5), MI_asynadapt_Mask(it,:), '-o')
    title(['Mask = ' num2str(testCons(it))])
    xlabel('Test contrast')
    ylabel('MI')
end
suptitle('Mask Pref')
