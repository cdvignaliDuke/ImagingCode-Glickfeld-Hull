edit('FlashingStim_dataSortedByCycle_combineDatasets.m')

dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_cycInd = cycV_ind{1};
A_cycInd = cycA_ind{1};
AV_cycInd = cycAV_ind{1};

%% find average response to first stim

dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_cycInd = cycV_ind{1};
A_cycInd = cycA_ind{1};
AV_cycInd = cycAV_ind{1};

V_cellAvgResp = squeeze(mean(dataTrialStart(35:41,:,V_cycInd),1))';
A_cellAvgResp = squeeze(mean(dataTrialStart(35:41,:,A_cycInd),1))';
AV_cellAvgResp = squeeze(mean(dataTrialStart(35:41,:,AV_cycInd),1))';

%% find max response for normalization - first stim

[maxResp_V, rInd_V] = max(V_cellAvgResp,[],1);
[maxResp_V, cInd_V] = max(maxResp_V,[],2);
rInd_V = rInd_V(1,cInd_V);

[maxResp_A, rInd_A] = max(A_cellAvgResp,[],1);
[maxResp_A, cInd_A] = max(maxResp_A,[],2);
rInd_A = rInd_A(1,cInd_A);

[maxResp_AV, rInd_AV] = max(AV_cellAvgResp,[],1);
[maxResp_AV, cInd_AV] = max(maxResp_AV,[],2);
rInd_AV = rInd_AV(1,cInd_AV);

maxResp = max([maxResp_V maxResp_A maxResp_AV]);
%% normalize responses - first stim

normCellResp_V = V_cellAvgResp/maxResp;
normCellResp_A = A_cellAvgResp/maxResp;
normCellResp_AV = AV_cellAvgResp/maxResp;

%% plot heatmap of responses across trials and accross cells - first stim
figure;
subplot(1,3,1)
imagesc(normCellResp_V(1:59,:));
hold on
title('visual trials');
xlabel('cells')
ylabel('trials')
colorbar('Location','SouthOutside')
caxis([0 1]);
subplot(1,3,2)
imagesc(normCellResp_A(1:59,:));
hold on
title('auditory trials');
xlabel('cells')
ylabel('trials')
colorbar('Location','SouthOutside')
caxis([0 1]);
subplot(1,3,3)
imagesc(normCellResp_AV(1:59,:));
hold on
title('vis+aud trials');
xlabel('cells')
ylabel('trials')
colorbar('Location','SouthOutside')
caxis([0 1]);



%% find average response to fourth stim

dataCycle4 = cycDataDFoverF_cmlvNoTarget{3};
V_cycInd4 = cycV_ind{3};
A_cycInd4 = cycA_ind{3};
AV_cycInd4 = cycAV_ind{3};

V_cellAvgResp4 = squeeze(mean(dataCycle4(68:73,:,V_cycInd4),1))';
A_cellAvgResp4 = squeeze(mean(dataCycle4(68:73,:,A_cycInd4),1))';
AV_cellAvgResp4 = squeeze(mean(dataCycle4(68:73,:,AV_cycInd4),1))';

%% find max response for normalization - first stim

[maxResp_V4, rInd_V4] = max(V_cellAvgResp4,[],1);
[maxResp_V4, cInd_V4] = max(maxResp_V4,[],2);
rInd_V4 = rInd_V4(1,cInd_V4);

[maxResp_A4, rInd_A4] = max(A_cellAvgResp4,[],1);
[maxResp_A4, cInd_A4] = max(maxResp_A4,[],2);
rInd_A4 = rInd_A4(1,cInd_A4);

[maxResp_AV4, rInd_AV4] = max(AV_cellAvgResp4,[],1);
[maxResp_AV4, cInd_AV4] = max(maxResp_AV4,[],2);
rInd_AV4 = rInd_AV4(1,cInd_AV4);

maxResp4 = max([maxResp_V4 maxResp_A4 maxResp_AV4]);
%% normalize responses - first stim

normCellResp_V4 = V_cellAvgResp4/maxResp4;
normCellResp_A4 = A_cellAvgResp4/maxResp4;
normCellResp_AV4 = AV_cellAvgResp4/maxResp4;

%% plot heatmap of responses across trials and accross cells - first stim
figure;
subplot(1,3,1)
imagesc(normCellResp_V4(1:42,:)');
hold on
title('visual trials');
xlabel('trials')
ylabel('cells')
colorbar('Location','SouthOutside')
caxis([0 1]);
subplot(1,3,2)
imagesc(normCellResp_A4(1:42,:)');
hold on
title('auditory trials');
xlabel('trials')
ylabel('cells')
colorbar('Location','SouthOutside')
caxis([0 1]);
subplot(1,3,3)
imagesc(normCellResp_AV4(1:42,:)');
hold on
title('vis+aud trials');
xlabel('trials')
ylabel('cells')
colorbar('Location','SouthOutside')
caxis([0 1]);



