edit('FlashingStim_dataSortedByCycle_combineDatasets.m')
%%
dirFolder = '005';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
cd(fileSave);
% load('oriTuningPreferences.mat')
load('TuningPreferences.mat')


%%
% cell types
dataTrialStart = cycDataDFoverF_cmlvNoTarget{4};
v_ind = cycV_ind{4};
% a_ind = cycA_ind{1};
a_ind = cycAV_ind{4};


preStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial =1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        preStimResp_V(itrial,icell) = mean(dataTrialStart(1:30,icell,v_ind(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial = 1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        baselineStimResp_V(itrial,icell) = mean(dataTrialStart(36:end,icell,v_ind(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);


cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsSelectZero_dir = intersect(dirSlctvCells,cellsPrefZero);
cellsSelectZero_ori = intersect(oriSlctvCells,cellsPrefZero);
cellsSelectZero = union(cellsSelectZero_dir,cellsSelectZero_ori);
cellsPrefRespZero = intersect(baselineStimRespIndex_V,cellsPrefZero);

cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
cellsSelectNinety_dir = intersect(dirSlctvCells,cellsPrefNinety);
cellsSelectNinety_ori = intersect(oriSlctvCells,cellsPrefNinety);
cellsSelectNinety = union(cellsSelectNinety_dir,cellsSelectNinety_ori);
cellsPrefRespNinety = intersect(baselineStimRespIndex_V,cellsPrefNinety);

nCells = size(cycDataDFoverF_cmlvNoTarget{7},2);
oriSlctvCellsAll = union(oriSlctvCells,dirSlctvCells);
notSlctvCells = setdiff([1:nCells],oriSlctvCellsAll);
notRespCells = setdiff([1:nCells],baselineStimRespIndex_V);

%%
% plot cells that are orientation or direction selective, prefer zero
% % degrees orientation, and significantly respond to baseline stim
% figure;
% for icyc = 1:length(cycles)
%     dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
%     V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
%     AV_cycInd = cycAV_ind{icyc};
%     V_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,A_cycInd),3),2);
%     AV_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,AV_cycInd),3),2);
%     subplot(3,3,icyc);
%     plot(V_avg,'g');
%     hold on
%     plot(A_avg,'r');
%     hold on
%     plot(AV_avg,'m');
%     hold on
%     vline(30,'k')
%     hold on
%     for i = 1:cycles(icyc)-1
%         L = (i*cycTime)+30;
%         vline(L,'k:');
%         hold on
%     end
%     vline((cycles(icyc)*cycTime+30),'c');
%     hold on
%     title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis; ' num2str(length(AV_cycInd)) ' aud+vis'])
%     hold on
% end

%% noise correlations 
% 
%find mean response to baseline flash for each cell
dataTrialStart = cycDataDFoverF_cmlvNoTarget{7};
V_cycInd = cycV_ind{7};
% A_cycInd = cycA_ind{7};
AV_cycInd = cycAV_ind{7};

cellGroup = 1:nCells;
cellGroupName = 'all'

% V_cellAvgResp = squeeze(mean(dataTrialStart(end-29:end,cellGroup,V_cycInd),1))';
% % A_cellAvgResp = squeeze(mean(dataTrialStart(:,cellsPrefZero,A_cycInd),1))';
% AV_cellAvgResp = squeeze(mean(dataTrialStart(end-29:end,cellGroup,AV_cycInd),1))';

%sort data by ori pref
[sortedOriPref, sortedOriPref_ind] = sort(oriPref_ind(cellGroup,:));
dataSelectCells = dataTrialStart(:,cellGroup,:);
dataSortedByOriPref = dataSelectCells(:,sortedOriPref_ind,:);
% [sortedRespDirPref, sortedRespOriPref_ind] = sort(oriPref_ind(baselineStimRespIndex_V));
% sortedRespIndex_V = baselineStimRespIndex_V(sortedRespOriPref_ind);

V_cellAvgResp = squeeze(mean(dataSortedByOriPref(end-29:end,:,V_cycInd),1))';
% A_cellAvgResp = squeeze(mean(dataSortedByOriPref(:,cellsPrefZero,A_cycInd),1))';
AV_cellAvgResp = squeeze(mean(dataSortedByOriPref(end-29:end,:,AV_cycInd),1))';

%find correlation coefficient (r) for cells
rSC_V = corrcoef(V_cellAvgResp);
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV = corrcoef(AV_cellAvgResp);
rSC_sub = rSC_AV-rSC_V;

rSC_V = tril(rSC_V,-1);
rSC_AV = tril(rSC_AV,-1);
rSC_sub = tril(rSC_sub,-1);

avgrSC_V = mean(mean(rSC_V(rSC_V  ~= 0),2),1);
avgrSC_AV = mean(mean(rSC_AV(rSC_AV  ~= 0),2),1);
avgrSC_sub = mean(mean(rSC_sub(rSC_sub  ~= 0),2),1);

rSC_sub_positive = rSC_sub >= 0;
rSC_sub_negative = rSC_sub < 0;
rSC_sub_percentpos = ((sum(sum(rSC_sub_positive)))/(size(rSC_sub,1).*size(rSC_sub,2)))*100;
rSC_sub_percentneg = ((sum(sum(rSC_sub_negative)))/(size(rSC_sub,1).*size(rSC_sub,2)))*100;

oriGroupBorders = zeros(size(unique(sortedOriPref)))';
runsum = 0;
for i = 1:length(oriGroupBorders)
    x = sum(sortedOriPref == i);
    oriGroupBorders(1,i) = x+runsum;
    runsum = runsum+x;
end

%plot noise correlation
figure;
colormap(brewermap([],'*RdBu'))
subplot(1,3,1);
imagesc(rSC_V)
hold on
title(['rSC_V avg corr = ' num2str(avgrSC_V)])
hold on
colorbar 
caxis([-1 1]);
axis('square')
hold on
for i = 1:length(oriGroupBorders)
    vline(oriGroupBorders(i),'k')
end
for i = 1:length(oriGroupBorders)
    hline(oriGroupBorders(i),'k')
end
set(gca,'XTick',oriGroupBorders)
set(gca,'XTickLabel',oriGroupBorders)
% subplot(1,3,2);
% imagesc(rSC_A)
% title('rSC_A')
% hold on
% colorbar
% caxis([-1 1]);
% axis('square')
subplot(1,3,2);
imagesc(rSC_AV)
title(['rSC_A_V avg corr = ' num2str(avgrSC_AV)])
hold on
colorbar
caxis([-1 1]);
axis('square')
hold on
for i = 1:length(oriGroupBorders)
    vline(oriGroupBorders(i),'k')
end
for i = 1:length(oriGroupBorders)
    hline(oriGroupBorders(i),'k')
end
set(gca,'XTick',oriGroupBorders)
set(gca,'XTickLabel',oriGroupBorders)
subplot(1,3,3);
imagesc(rSC_sub)
hold on
title(['AV - V subtraction, avg corr = ' num2str(avgrSC_sub)])
hold on
colorbar
caxis([-1 1]);
axis('square')
hold on
for i = 1:length(oriGroupBorders)
    vline(oriGroupBorders(i),'k')
end
for i = 1:length(oriGroupBorders)
    hline(oriGroupBorders(i),'k')
end
set(gca,'XTick',oriGroupBorders)
set(gca,'XTickLabel',oriGroupBorders)
ha = axes('Position',[0 0 1 1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');
text(0.5, 1,['\bf ' mouse ' ' date '; ' cellGroupName ' cells'],'HorizontalAlignment','center','VerticalAlignment', 'top')

%%
%bin cells by orientation preference
% 
%find mean response to baseline flash for each cell
dataTrialStart = cycDataDFoverF_cmlvNoTarget{9};
V_cycInd = cycV_ind{9};
% A_cycInd = cycA_ind{7};
AV_cycInd = cycAV_ind{9};

cellGroup = baselineStimRespIndex_V;

% V_cellAvgResp = squeeze(mean(dataTrialStart(end-29:end,cellGroup,V_cycInd),1))';
% % A_cellAvgResp = squeeze(mean(dataTrialStart(:,cellsPrefZero,A_cycInd),1))';
% AV_cellAvgResp = squeeze(mean(dataTrialStart(end-29:end,cellGroup,AV_cycInd),1))';

%sort data by ori pref
[sortedOriPref, sortedOriPref_ind] = sort(oriPref_ind(cellGroup,:));
dataSelectCells = dataTrialStart(:,cellGroup,:);
dataSortedByOriPref = dataSelectCells(:,sortedOriPref_ind,:);
% [sortedRespDirPref, sortedRespOriPref_ind] = sort(oriPref_ind(baselineStimRespIndex_V));
% sortedRespIndex_V = baselineStimRespIndex_V(sortedRespOriPref_ind);

V_cellAvgResp = squeeze(mean(dataSortedByOriPref(end-29:end,:,V_cycInd),1))';
% A_cellAvgResp = squeeze(mean(dataSortedByOriPref(:,cellsPrefZero,A_cycInd),1))';
AV_cellAvgResp = squeeze(mean(dataSortedByOriPref(end-29:end,:,AV_cycInd),1))';

%find correlation coefficient (r) for cells
rSC_V = corrcoef(V_cellAvgResp);
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV = corrcoef(AV_cellAvgResp);
rSC_sub = rSC_AV-rSC_V;

rSC_V = tril(rSC_V,-1);
rSC_AV = tril(rSC_AV,-1);
rSC_sub = tril(rSC_sub,-1);

avgrSC_V = mean(mean(rSC_V(rSC_V  ~= 0),2),1);
avgrSC_AV = mean(mean(rSC_AV(rSC_AV  ~= 0),2),1);
avgrSC_sub = mean(mean(rSC_sub(rSC_sub  ~= 0),2),1);

rSC_sub_positive = rSC_sub >= 0;
rSC_sub_negative = rSC_sub < 0;
rSC_sub_percentpos = ((sum(sum(rSC_sub_positive)))/(size(rSC_sub,1).*size(rSC_sub,2)))*100;
rSC_sub_percentneg = ((sum(sum(rSC_sub_negative)))/(size(rSC_sub,1).*size(rSC_sub,2)))*100;

oriGroupBorders = zeros(size(unique(sortedOriPref)))';
runsum = 0;
for i = 1:length(oriGroupBorders)
    x = sum(sortedOriPref == i);
    oriGroupBorders(1,i) = x+runsum;
    runsum = runsum+x;
end


for i = 1:4
    ind = find(sortedOriPref == i);
    logicInd = zeros(1,length(oriPref_ind));
    logicInd(:,ind) = 1;
    cellOriPref_ind(i,:) = logicInd;
end
cellOriPref_ind = logical(cellOriPref_ind);

for i = 1:4
    V_respPrefAvg(i) = mean(V_avg(cellOriPref_ind(i,:)));
    AV_respPrefAvg(i) = mean(AV_avg(cellOriPref_ind(i,:)));
end

rSC_Vvector = rSC_V(:);
rSC_Vvector = rSC_Vvector(rSC_Vvector ~= 0);
rSC_AVvector = rSC_AV(:);
rSC_AVvector = rSC_AVvector(rSC_AVvector ~= 0);
rSC_subVector = rSC_sub(:);
rSC_subVector = rSC_subVector(rSC_subVector ~= 0);
for i = 1:4
    rSC_VrespPrefAvg(i) = mean(rSC_Vvector(cellOriPref_ind(i,:)));
    rSC_AVrespPrefAvg(i) = mean(rSC_AVvector(cellOriPref_ind(i,:)));
end

for i = 1:4
    errbar_rSCV_ori(i) = (std(rSC_Vavg(cellOriPref_ind(i,:))))/(sqrt(size(rSC_Vavg(cellOriPref_ind(i,:)),1)));
    errbar_rSCAV_ori(i) = (std(rSC_AVavg(cellOriPref_ind(i,:))))/(sqrt(size(rSC_AVavg(cellOriPref_ind(i,:)),1)));
end

figure;
errorbar(rSC_VrespPrefAvg,errbar_rSCV_ori,'g');
hold on
errorbar(rSC_AVrespPrefAvg,errbar_rSCAV_ori,'m')
hold on
xlim([0 5])
ylim([0.05 0.15])
xlabel('Orientation Preference')
ylabel('noise corr')
set(gca,'xTick', [1 2 3 4])
set(gca,'xTickLabel', {'0' '45' '90' '135'})


%%
% other plots
V_cellAvgResp = squeeze(mean(dataTrialStart(:,:,V_cycInd),1))';
% A_cellAvgResp = squeeze(mean(dataTrialStart(:,cellsSelectZero,A_cycInd),1))';
AV_cellAvgResp = squeeze(mean(dataTrialStart(:,:,AV_cycInd),1))';

%find correlation coefficient (r) for cells
rSC_V = corrcoef(V_cellAvgResp);
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV = corrcoef(AV_cellAvgResp);

rSC_sub = rSC_AV-rSC_V;
rSC_sub_positive = rSC_sub >= 0;
rSC_sub_negative = rSC_sub < 0;
rSC_sub_percentpos = ((sum(sum(rSC_sub_positive)))/(size(rSC_sub,1).*size(rSC_sub,2)))*100;
rSC_sub_percentneg = ((sum(sum(rSC_sub_negative)))/(size(rSC_sub,1).*size(rSC_sub,2)))*100;

avgrSC_V = mean(mean(rSC_V,2),1);
avgrSC_AV = mean(mean(rSC_AV,2),1);
avgrSC_sub = mean(mean(rSC_sub,2),1);

% bin cells by average firing rate
V_avg = squeeze(mean(V_cellAvgResp,1));
AV_avg = squeeze(mean(AV_cellAvgResp,1));

minFR = min([V_avg AV_avg]);
maxFR = max([V_avg AV_avg]);
binedges = linspace(minFR,maxFR,10);
binedges = [binedges(1:end-1) 1];

for i = 1:length(binedges(1:end-1))
    V_binInd(i,:) = V_avg >= binedges(i) & V_avg < binedges(i+1);
    AV_binInd(i,:) = AV_avg >= binedges(i) & AV_avg < binedges(i+1);
end

cellsPerVbin = sum(V_binInd,2);
cellsPerAVbin = sum(AV_binInd,2);

for i = 1:length(binedges(1:end-1))
    V_binAvg(i) = mean(V_avg(V_binInd(i,:)));
    AV_binAvg(i) = mean(AV_avg(AV_binInd(i,:)));
end

rSC_Vavg = mean(rSC_V,2);
rSC_AVavg = mean(rSC_AV,2);
for i = 1:length(binedges(1:end-1))
    rSC_VbinAvg(i) = mean(rSC_Vavg(V_binInd(i,:)));
    rSC_AVbinAvg(i) = mean(rSC_AVavg(AV_binInd(i,:)));
end

for i = 1:length(binedges(1:end-1))
    errbar_rSCV(i) = (std(rSC_Vavg(V_binInd(i,:))))/(sqrt(size(rSC_Vavg(V_binInd(i,:)),1)));
    errbar_rSCAV(i) = (std(rSC_AVavg(AV_binInd(i,:))))/(sqrt(size(rSC_AVavg(AV_binInd(i,:)),1)));
end

figure;
errorbar(V_binAvg,rSC_VbinAvg,errbar_rSCV,'g');
hold on
errorbar(AV_binAvg,rSC_AVbinAvg,errbar_rSCAV,'m')
xlabel('dF/F')
ylabel('noise corr')


%% noise correlations - all cells
% 
%find mean response to baseline flash for each cell
dataTrialStart = cycDataDFoverF_cmlvNoTarget{6};
V_cycInd = cycV_ind{6};
% A_cycInd = cycA_ind{1};
AV_cycInd = cycAV_ind{6};


V_cellAvgResp = squeeze(mean(dataTrialStart(:,:,V_cycInd),1))';
% A_cellAvgResp = squeeze(mean(dataTrialStart(:,cellsSelectZero,A_cycInd),1))';
AV_cellAvgResp = squeeze(mean(dataTrialStart(:,:,AV_cycInd),1))';

%find correlation coefficient (r) for cells
rSC_V = corrcoef(V_cellAvgResp);
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV = corrcoef(AV_cellAvgResp);

rSC_sub = rSC_AV-rSC_V;
rSC_sub_positive = rSC_sub >= 0;
rSC_sub_negative = rSC_sub < 0;
rSC_sub_percentpos = ((sum(sum(rSC_sub_positive)))/(size(rSC_sub,1).*size(rSC_sub,2)))*100;
rSC_sub_percentneg = ((sum(sum(rSC_sub_negative)))/(size(rSC_sub,1).*size(rSC_sub,2)))*100;

avgrSC_V = mean(mean(rSC_V,2),1);
avgrSC_AV = mean(mean(rSC_AV,2),1);
avgrSC_sub = mean(mean(rSC_sub,2),1);

% bin cells by average firing rate
V_avg = squeeze(mean(V_cellAvgResp,1));
AV_avg = squeeze(mean(AV_cellAvgResp,1));

minFR = min([V_avg AV_avg]);
maxFR = max([V_avg AV_avg]);
binedges = linspace(minFR,maxFR,10);
binedges = [binedges(1:end-1) 1];

for i = 1:length(binedges(1:end-1))
    V_binInd(i,:) = V_avg >= binedges(i) & V_avg < binedges(i+1);
    AV_binInd(i,:) = AV_avg >= binedges(i) & AV_avg < binedges(i+1);
end

cellsPerVbin = sum(V_binInd,2);
cellsPerAVbin = sum(AV_binInd,2);

for i = 1:length(binedges(1:end-1))
    V_binAvg(i) = mean(V_avg(V_binInd(i,:)));
    AV_binAvg(i) = mean(AV_avg(AV_binInd(i,:)));
end

rSC_Vavg = mean(rSC_V,2);
rSC_AVavg = mean(rSC_AV,2);
for i = 1:length(binedges(1:end-1))
    rSC_VbinAvg(i) = mean(rSC_Vavg(V_binInd(i,:)));
    rSC_AVbinAvg(i) = mean(rSC_AVavg(AV_binInd(i,:)));
end

figure;
plot(V_binAvg,rSC_VbinAvg,'g');
hold on
plot(AV_binAvg,rSC_AVbinAvg,'m')
xlabel('dF/F')
ylabel('noise corr')

%bin cells by orientation preference
for i = 1:4
    ind = find(dirPref_ind == i | dirPref_ind == i+4);
    logicInd = zeros(1,length(dirPref_ind));
    logicInd(:,ind) = 1;
    cellOriPref_ind(i,:) = logicInd;
end
cellOriPref_ind = logical(cellOriPref_ind);

for i = 1:4
    V_respPrefAvg(i) = mean(V_avg(cellOriPref_ind(i,:)));
    AV_respPrefAvg(i) = mean(AV_avg(cellOriPref_ind(i,:)));
end

rSC_Vavg = mean(rSC_V,2);
rSC_AVavg = mean(rSC_AV,2);
for i = 1:4
    rSC_VrespPrefAvg(i) = mean(rSC_Vavg(cellOriPref_ind(i,:)));
    rSC_AVrespPrefAvg(i) = mean(rSC_AVavg(cellOriPref_ind(i,:)));
end

figure;
plot(rSC_VrespPrefAvg,'g');
hold on
plot(rSC_AVrespPrefAvg,'m')
xlim([0 5])
ylim([0.05 0.15])
xlabel('Orientation Preference')
ylabel('noise corr')
set(gca,'xTick', [1 2 3 4])
set(gca,'xTickLabel', {'0' '45' '90' '135'})



%plot signal correlation
figure;
colormap(brewermap([],'*RdBu'))
subplot(1,3,1);
imagesc(rSC_V)
hold on
title('rSC_V')
hold on
colorbar 
caxis([-1 1]);
axis('square')
hold on
subplot(1,3,2);
imagesc(rSC_AV)
title('rSC_A_V')
hold on
colorbar
caxis([-1 1]);
axis('square')
subplot(1,3,3);
imagesc(rSC_sub)
hold on
title(['subtraction, ' num2str(rSC_sub_percentpos) '% pos'])
hold on
colorbar
caxis([-1 1]);
axis('square')

%% plot correlations for all selectivity types
[sortedFR, sortedFR_ind] = sort(squeeze(mean(V_cellAvgResp,1)));
dataSortedByFR = dataTrialStart(:,cellsPrefZero,:);
dataSortedByFR = dataSortedByFR(:,sortedFR_ind,:);

V_cellAvgResp = squeeze(mean(dataSortedByFR(end-19:end,:,V_cycInd),1))';
% A_cellAvgResp = squeeze(mean(dataSortedByFR(35:41,cellsSelectZero,A_cycInd),1))';
AV_cellAvgResp = squeeze(mean(dataSortedByFR(end-19:end,:,AV_cycInd),1))';

rSC_V = corrcoef(V_cellAvgResp);
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV = corrcoef(AV_cellAvgResp);


%plot noise corr x FR
figure;
plot(squeeze(mean(V_cellAvgResp)),rSC_V,'g')
hold on
plot(squeeze(mean(AV_cellAvgResp)),rSC_AV,'g')

%plot signal correlation
figure;
subplot(1,3,1);
imagesc(rSC_V)
hold on
title('rSC_V')
hold on
colorbar
caxis([0 1]);
hold on
subplot(1,3,2);
imagesc(rSC_AV)
title('rSC_A_V')
hold on
colorbar
caxis([0 1]);
subplot(1,3,3);
imagesc(rSC_AV-rSC_V)
hold on
title('sub')
hold on
colorbar
caxis([0 1]);

%% rSC - all cells
dataTrialStart = cycDataDFoverF_cmlvNoTarget{7};
V_cycInd = cycV_ind{7};
% A_cycInd = cycA_ind{7};
AV_cycInd = cycAV_ind{7};
[sortedDirPref, sortedDirPref_ind] = sort(dirPref_ind);
dataSortedByDirPref = dataTrialStart(:,sortedDirPref_ind,:);
[sortedRespDirPref, sortedRespDirPref_ind] = sort(dirPref_ind(baselineStimRespIndex_V));
sortedRespIndex_V = baselineStimRespIndex_V(sortedRespDirPref_ind);

V_cellAvgResp_all = squeeze(mean(dataSortedByDirPref(:,sortedRespIndex_V,V_cycInd),1))';
% A_cellAvgResp_all = squeeze(mean(dataSortedByDirPref(:,sortedRespIndex_V,A_cycInd),1))';
AV_cellAvgResp_all = squeeze(mean(dataSortedByDirPref(:,sortedRespIndex_V,AV_cycInd),1))';
% V_cellAvgResp_all = squeeze(mean(dataTrialStart(:,baselineStimRespIndex_V,V_cycInd),1))';
% % A_cellAvgResp_all = squeeze(mean(dataTrialStart(:,baselineStimRespIndex_V,A_cycInd),1))';
% AV_cellAvgResp_all = squeeze(mean(dataTrialStart(:,baselineStimRespIndex_V,AV_cycInd),1))';

rSC_V_all = corrcoef(V_cellAvgResp_all);
% rSC_A_all = corrcoef(A_cellAvgResp_all);
rSC_AV_all = corrcoef(AV_cellAvgResp_all);
rSC_sub_all = rSC_AV_all - rSC_V_all;

avgrSC_V_all = mean(mean(rSC_V_all,2),1);
avgrSC_AV_all = mean(mean(rSC_AV_all,2),1);
avgrSC_sub_all = mean(mean(rSC_sub_all,2),1);


%plot signal correlation
figure;
colormap(brewermap([],'*RdBu'))
subplot(1,3,1);
imagesc(rSC_V_all)
hold on
title(['rSC_V avg corr = ' num2str(avgrSC_V_all)])
hold on
colorbar 
caxis([-1 1]);
axis('square')
hold on
subplot(1,3,2);
imagesc(rSC_AV_all)
title(['rSC_A_V avg corr = ' num2str(avgrSC_AV_all)])
hold on
colorbar
caxis([-1 1]);
axis('square')
subplot(1,3,3);
imagesc(rSC_sub_all)
hold on
title(['AV - V subtraction, avg corr = ' num2str(avgrSC_sub_all)])
hold on
colorbar
caxis([-1 1]);
axis('square')

figure;
subplot(1,3,1);
imagesc(rSC_V_all)
hold on
title('rSC_V')
hold on
colorbar('Location','SouthOutside');
caxis([0 1]);
hold on

subplot(1,3,2);
imagesc(rSC_AV_all)
title('rSC_A_V')
hold on
colorbar('Location','SouthOutside')
caxis([0 1]);
subplot(1,3,3);
imagesc(rSC_AV_all-rSC_V_all)
hold on
title('sub')
hold on
colorbar('Location','SouthOutside')
caxis([0 1]);

rSC_sub_all = rSC_AV_all-rSC_V_all;
rSC_sub_positive = rSC_sub_all >= 0;
rSC_sub_negative = rSC_sub_all < 0;
rSC_sub_percentpos = (sum(sum(rSC_sub_positive)))/(size(rSC_sub_all,1).*size(rSC_sub_all,2));
rSC_sub_percentneg = (sum(sum(rSC_sub_negative)))/(size(rSC_sub_all,1).*size(rSC_sub_all,2));

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
