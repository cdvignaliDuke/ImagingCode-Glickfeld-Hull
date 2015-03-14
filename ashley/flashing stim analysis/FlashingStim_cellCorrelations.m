edit('FlashingStim_dataSortedByCycle_combineDatasets.m')

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('cellSelectivityIndices.mat')

cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsSelectZero_dir = intersect(dirSlctvCells,cellsPrefZero);
cellsSelectZero_ori = intersect(oriSlctvCells,cellsPrefZero);
cellsSelectZero = union(cellsSelectZero_dir,cellsSelectZero_ori);
respCellsSelectZero = intersect(baselineStimRespIndex_V,cellsSelectZero);

%%
% plot cells that are orientation or direction selective, prefer zero
% degrees orientation, and significantly respond to baseline stim
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,V_cycInd),3),2);
    A_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(A_avg,'r');
    hold on
    plot(AV_avg,'m');
    hold on
    vline(30,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis; ' num2str(length(AV_cycInd)) ' aud+vis'])
    hold on
end

%% "spike count" correlations - cells selective for zero
% 
%find mean response to baseline flash for each cell
plotRows = ceil(size(respCellsSelectZero,1)/4);
dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_cycInd = cycV_ind{1};
A_cycInd = cycA_ind{1};
AV_cycInd = cycAV_ind{1};


V_cellAvgResp = squeeze(mean(dataTrialStart(35:41,respCellsSelectZero,V_cycInd),1))';
A_cellAvgResp = squeeze(mean(dataTrialStart(35:41,respCellsSelectZero,A_cycInd),1))';
AV_cellAvgResp = squeeze(mean(dataTrialStart(35:41,respCellsSelectZero,AV_cycInd),1))';

%find correlation coefficient (r) for cells
rSC_V = corrcoef(V_cellAvgResp);
rSC_A = corrcoef(A_cellAvgResp);
rSC_AV = corrcoef(AV_cellAvgResp);

%plot signal correlation
figure;
subplot(1,3,1);
imagesc(rSC_V)
hold on
title('rSC_V')
hold on
colorbar
hold on
subplot(1,3,2);
imagesc(rSC_A)
hold on
title('rSC_A')
hold on
colorbar
hold on
subplot(1,3,3);
imagesc(rSC_AV)
title('rSC_A_V')
hold on
colorbar

%% rSC - all cells
[sortedDirPref, sortedDirPref_ind] = sort(dirPref_ind);
dataSortedByDirPref = dataTrialStart(:,sortedDirPref_ind,:);

V_cellAvgResp_all = squeeze(mean(dataSortedByDirPref(35:41,:,V_cycInd),1))';
A_cellAvgResp_all = squeeze(mean(dataSortedByDirPref(35:41,:,A_cycInd),1))';
AV_cellAvgResp_all = squeeze(mean(dataSortedByDirPref(35:41,:,AV_cycInd),1))';

rSC_V_all = corrcoef(V_cellAvgResp_all);
rSC_A_all = corrcoef(A_cellAvgResp_all);
rSC_AV_all = corrcoef(AV_cellAvgResp_all);

figure;
subplot(1,3,1);
imagesc(rSC_V_all)
hold on
title('rSC_V')
hold on
colorbar('Location','SouthOutside');
caxis([-.3 1]);
hold on
subplot(1,3,2);
imagesc(rSC_A_all)
hold on
title('rSC_A')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);
hold on
subplot(1,3,3);
imagesc(rSC_AV_all)
title('rSC_A_V')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);



dirs = unique(dirPref_ind);
count = 0;
for i = 1:size(dirs,1);
    numDirPref(i) = sum(dirPref_ind == dirs(i))+count;
    count = count+sum(dirPref_ind == dirs(i));
end

%% correlations at 4 cycles
[sortedDirPref, sortedDirPref_ind] = sort(dirPref_ind);
data4Cyc = cycDataDFoverF_cmlvNoTarget{3};
dataSortedByDirPref = data4Cyc(:,sortedDirPref_ind,:);

V_cyc4Ind = cycV_ind{3};
A_cyc4Ind = cycA_ind{3};
AV_cyc4Ind = cycAV_ind{3};

% V_cellAvgResp = squeeze(mean(data4Cyc(68:74,respCellsSelectZero,V_cyc4Ind),1))';
% % A_cellAvgResp = squeeze(mean(data4Cyc(68:74,respCellsSelectZero,A_cyc4Ind),1))';
% AV_cellAvgResp = squeeze(mean(data4Cyc(68:74,respCellsSelectZero,AV_cyc4Ind),1))';
% 
% %find correlation coefficient (r) for cells
% rSC_V = corrcoef(V_cellAvgResp);
% % rSC_A = corrcoef(A_cellAvgResp);
% rSC_AV = corrcoef(AV_cellAvgResp);
% 
% %plot signal correlation
% figure;
% subplot(1,3,1);
% imagesc(rSC_V)
% hold on
% title('rSC_V')
% hold on
% colorbar
% hold on
% % subplot(1,3,2);
% % imagesc(rSC_A)
% % hold on
% % title('rSC_A')
% % hold on
% % colorbar
% % hold on
% subplot(1,3,3);
% imagesc(rSC_AV)
% title('rSC_A_V')
% hold on
% colorbar

%all cells

V_cellAvgResp_all_cyc4 = squeeze(mean(dataSortedByDirPref(68:74,:,V_cyc4Ind),1))';
A_cellAvgResp_all = squeeze(mean(dataSortedByDirPref(68:74,:,A_cyc4Ind),1))';
AV_cellAvgResp_all_cyc4 = squeeze(mean(dataSortedByDirPref(68:74,:,AV_cyc4Ind),1))';

rSC_V_all_cyc4 = corrcoef(V_cellAvgResp_all_cyc4);
rSC_A_all = corrcoef(A_cellAvgResp_all);
rSC_AV_all_cyc4 = corrcoef(AV_cellAvgResp_all_cyc4);

figure;
subplot(1,3,1);
imagesc(rSC_V_all_cyc4)
hold on
title('rSC_V')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);
hold on
subplot(1,3,2);
imagesc(rSC_A_all)
hold on
title('rSC_A')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);
hold on
subplot(1,3,3);
imagesc(rSC_AV_all_cyc4)
title('rSC_A_V')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);

% figure;
% subplot(1,4,1)
% imagesc(rSC_V_all)
% hold on
% title('rSC_V - cyc1')
% hold on
% colorbar
% hold on
% subplot(1,4,2);
% imagesc(rSC_AV_all)
% title('rSC_A_V - cyc1')
% hold on
% colorbar
% hold on
% subplot(1,4,3);
% imagesc(rSC_V_all_cyc4)
% hold on
% title('rSC_V - cyc4')
% hold on
% colorbar
% hold on
% subplot(1,4,4);
% imagesc(rSC_AV_all_cyc4)
% title('rSC_A_V - cyc4')
% hold on
% colorbar

%% subtract rSignal_V from rSignal_AV to see difference in correlations

rSC_subV = rSC_AV_all - rSC_V_all;
rSC_subAV = rSC_V_all - rSC_AV_all;
rSC_V_subA = rSC_V_all - rSC_A_all;
rSC_AV_subA = rSC_AV_all - rSC_A_all;

%plot signal correlation
figure;
subplot(1,5,1);
imagesc(rSC_V_all)
hold on
title('rSC_V')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);
hold on
subplot(1,5,2);
imagesc(rSC_A_all)
hold on
title('rSC_A')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);
hold on
subplot(1,5,3);
imagesc(rSC_AV_all)
title('rSC_A_V')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);
hold on
subplot(1,5,4);
imagesc(rSC_AV_subA)
title('rSC-subV')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);
hold on
subplot(1,5,5);
imagesc(rSC_AV_subAandV)
title('rSC-subAV')
hold on
colorbar('Location','SouthOutside')
caxis([-.3 1]);
