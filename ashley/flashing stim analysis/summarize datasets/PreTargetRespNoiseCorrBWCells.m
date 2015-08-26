% run(FlashingStim_dataSortedByCycle_combineDatasets.m)
%%
% oriselectivecells, driven by baseline stim, "driven
% cells", non-selective cells, target driven cells
%%
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' date '\PreTargetNoiseCorr'];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' date]);
        mkdir('PreTargetNoiseCorr')
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(date,'PreTargetNoiseCorr')
    end
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
%% find sets of cells
DirFolder = '006';
run('cellSets.m')
%%
cells = oriSlctvCellsAll;
cellgroupname = 'ori selective cells';
figBaseName = 'noiseCorr_orislctvcells_success';
%% set data
CYC = find(cycles == 8);
data = cycDataDFoverF_cmlvNoTarget{CYC};
V_cycInd = intersect(cycV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));
AV_cycInd = intersect(cycAV_ind{CYC},find(strcmp(cycTrialOutcome{CYC},'success')));

%sort data by ori pref
[sortedOriPref, sortedOriPref_ind] = sort(oriPref_ind(cells,:));
dataSelectCells = data(:,cells,:);
dataSortedByOriPref = dataSelectCells(:,sortedOriPref_ind,:);

%% noise correlations last 3cyles of long trials

V_cellAvgRespEND = squeeze(mean(dataSortedByOriPref(end-(cycTime*3):end,:,V_cycInd),1))';
AV_cellAvgRespEND = squeeze(mean(dataSortedByOriPref(end-(cycTime*3):end,:,AV_cycInd),1))';

%find correlation coefficient (r) for cells
rSC_V_END = corrcoef(V_cellAvgRespEND);
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV_END = corrcoef(AV_cellAvgRespEND);
rSC_sub_END = rSC_AV_END-rSC_V_END;

rSC_V_END = tril(rSC_V_END,-1);
rSC_AV_END = tril(rSC_AV_END,-1);
rSC_sub_END = tril(rSC_sub_END,-1);

avgrSC_V_END = mean(mean(rSC_V_END(rSC_V_END  ~= 0),2),1);
avgrSC_AV_END = mean(mean(rSC_AV_END(rSC_AV_END  ~= 0),2),1);
avgrSC_sub_END = mean(mean(rSC_sub_END(rSC_sub_END  ~= 0),2),1);

%% noise correlations last 3cyles of long trials
V_cellAvgRespSTART = squeeze(mean(dataSortedByOriPref(31:30+(cycTime*3),:,V_cycInd),1))';
AV_cellAvgRespSTART = squeeze(mean(dataSortedByOriPref(31:30+(cycTime*3),:,AV_cycInd),1))';

%find correlation coefficient (r) for cells
rSC_V_START = corrcoef(V_cellAvgRespSTART);
% rSC_A = corrcoef(A_cellAvgResp);
rSC_AV_START = corrcoef(AV_cellAvgRespSTART);
rSC_sub_START = rSC_AV_START-rSC_V_START;

rSC_V_START = tril(rSC_V_START,-1);
rSC_AV_START = tril(rSC_AV_START,-1);
rSC_sub_START = tril(rSC_sub_START,-1);

avgrSC_V_START = mean(mean(rSC_V_START(rSC_V_START  ~= 0),2),1);
avgrSC_AV_START = mean(mean(rSC_AV_START(rSC_AV_START  ~= 0),2),1);
avgrSC_sub_START = mean(mean(rSC_sub_START(rSC_sub_START  ~= 0),2),1);

%% plot noise correlations first 3cycles vs last 3cyles of long trials
oriGroupBorders = zeros(size(unique(sortedOriPref)))';
runsum = 0;
for i = 1:length(oriGroupBorders)
    x = sum(sortedOriPref == i);
    oriGroupBorders(1,i) = x+runsum;
    runsum = runsum+x;
end

%plot noise correlation
figure;
suptitle({[mouse '; ' date '; ' cellgroupname 'noise corr'];['top = end of trial, bottom = beginning'];[num2str(length(V_cycInd)) ' vis, ' num2str(length(AV_cycInd)) 'aud; ' num2str(length(cells)) ' cells']})
colormap(brewermap([],'*RdBu'))
subplot(2,3,1);
imagesc(rSC_V_END)
hold on
title(['rSC_V avg corr = ' num2str(avgrSC_V_END)])
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
subplot(2,3,2);
imagesc(rSC_AV_END)
title(['rSC_A_V avg corr = ' num2str(avgrSC_AV_END)])
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
subplot(2,3,3);
imagesc(rSC_sub_END)
hold on
title(['AV - V subtraction, avg corr = ' num2str(avgrSC_sub_END)])
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

colormap(brewermap([],'*RdBu'))
subplot(2,3,4);
imagesc(rSC_V_START)
hold on
title(['rSC_V avg corr = ' num2str(avgrSC_V_START)])
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
subplot(2,3,5);
imagesc(rSC_AV_START)
title(['rSC_A_V avg corr = ' num2str(avgrSC_AV_START)])
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
subplot(2,3,6);
imagesc(rSC_sub_START)
hold on
title(['AV - V subtraction, avg corr = ' num2str(avgrSC_sub_START)])
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

print([fnout ['\' figBaseName '_heatmap' '.pdf']], '-dpdf')

%% plot avg noise corr by ori pref
avgrSC_V_START_byOri = zeros(4,1);
avgrSC_V_END_byOri = zeros(4,1);
errrSC_V_START_byOri = zeros(4,1);
errrSC_V_END_byOri = zeros(4,1);
avgrSC_AV_START_byOri = zeros(4,1);
avgrSC_AV_END_byOri = zeros(4,1);
errrSC_AV_START_byOri = zeros(4,1);
errrSC_AV_END_byOri = zeros(4,1);
for i = 1:4
    ind = find(sortedOriPref == i);
    tempdataVSTART = rSC_V_START(ind,ind);
    tempdataVEND = rSC_V_END(ind,ind);
    avgrSC_V_START_byOri(i,:) = mean(tempdataVSTART(tempdataVSTART ~= 0),1);
    avgrSC_V_END_byOri(i,:) = mean(tempdataVEND(tempdataVEND ~= 0),1);
    errrSC_V_START_byOri(i,:) = (std(tempdataVSTART(tempdataVSTART ~= 0),0,1))/(sqrt(length(tempdataVSTART(tempdataVSTART ~= 0))));
    errrSC_V_END_byOri(i,:) = (std(tempdataVEND(tempdataVEND ~= 0),0,1))/(sqrt(length(tempdataVEND(tempdataVEND ~= 0))));
    
    tempdataAVSTART = rSC_AV_START(ind,ind);
    tempdataAVEND = rSC_AV_END(ind,ind);
    avgrSC_AV_START_byOri(i,:) = mean(tempdataAVSTART(tempdataAVSTART ~= 0),1);
    avgrSC_AV_END_byOri(i,:) = mean(tempdataAVEND(tempdataAVEND ~= 0),1);
    errrSC_AV_START_byOri(i,:) = (std(tempdataAVSTART(tempdataAVSTART ~= 0),0,1))/(sqrt(length(tempdataAVSTART(tempdataAVSTART ~= 0))));
    errrSC_AV_END_byOri(i,:) = (std(tempdataAVEND(tempdataAVEND ~= 0),0,1))/(sqrt(length(tempdataAVEND(tempdataAVEND ~= 0))));
end

figure;
suptitle({[mouse '; ' date '; ' cellgroupname 'noise corr'];['avg noise corr by ori pref'];[num2str(length(V_cycInd)) ' vis, ' num2str(length(AV_cycInd)) 'aud']})
subplot(1,2,1)
errorbar(avgrSC_V_START_byOri,errrSC_V_START_byOri,'g')
hold on
errorbar(avgrSC_AV_START_byOri,errrSC_AV_START_byOri,'k')
% xlim([0 5])
ylim([-0.3 0.3])
xlabel('Orientation Preference')
ylabel('noise corr')
title('beginning of trial')
set(gca,'xTick', [1 2 3 4])
set(gca,'xTickLabel', {'0' '45' '90' '135'})

subplot(1,2,2)
errorbar(avgrSC_V_END_byOri,errrSC_V_END_byOri,'g')
hold on
errorbar(avgrSC_AV_END_byOri,errrSC_AV_END_byOri,'k')
% xlim([0 5])
ylim([-0.3 0.3])
xlabel('Orientation Preference')
ylabel('noise corr')
title('end of trial')
set(gca,'xTick', [1 2 3 4])
set(gca,'xTickLabel', {'0' '45' '90' '135'})
legend({'vis'; 'vis+aud'},'Location','SouthWest')

print([fnout ['\' figBaseName '_avgbyoripref' '.pdf']], '-dpdf')