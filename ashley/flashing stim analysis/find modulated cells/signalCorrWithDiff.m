run('FlashingStim_dataSortedByCycle_combineDatasets.m')
%%
dirFolder = '006';
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

for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));   
    cycTrialOutcome{icyc} = trialOutcome(ind);
    cycDirectionDeg{icyc} = DirectionDeg(ind);
%     cycCatchDirectionDeg{icyc} = catchDirectionDeg(ind);
%     cycCatchTrialOutcome{icyc} = catchTrialOutcome(ind);
%     cycCatchCycle{icyc} = catchCycle(ind);
end
%%
dataavgtrials = squeeze(mean(dataDFoverF(:,:,V_cycInd),3));
drivencells = find(any(dataavgtrials >0.05,1));
%%
% cell1 = 2;
% tempdata1 = dataDFoverF;
% 
% % for i = 1:length(drivencells)
% figure;
% % 
% %     subplot(5,5,i)
% plot(tsmovavg(squeeze(tempdata1(:,(drivencells(cell1)),V_cycInd)),'s',3,1),'g')
% alpha(0.25)
% hold on
% plot(tsmovavg(squeeze(tempdata1(:,(drivencells(cell1)),AV_cycInd)),'s',3,1),'k')
% % hold on
% vline(30,'k') 
%     for iL = 1:cycles(icyc)-1
%         L = (iL*cycTime)+31;
%         vline(L,'k:');
%         hold on
%     end
% hold on
% plot(mean(tsmovavg(squeeze(tempdata1(:,(drivencells(cell1)),V_cycInd)),'s',3,1),2),'g','LineWidth',3)
% hold on
% plot(mean(tsmovavg(squeeze(tempdata1(:,(drivencells(cell1)),AV_cycInd)),'s',3,1),2),'k','LineWidth',3)
% % end

%%
V_cycInd = cycV_ind{end-1};
AV_cycInd = cycAV_ind{end-1};
datadiff = zeros(size(cycDataDFoverF_cmlvNoTarget{end-1}));
for icell = 1:size(cycDataDFoverF_cmlvNoTarget{end-1},2)
for itrial = 1:size(cycDataDFoverF_cmlvNoTarget{end-1},3)
    datadiff(2:end,icell,itrial) = diff(squeeze(cycDataDFoverF_cmlvNoTarget{end-1}(:,icell,itrial)));
end
end
datadiff(datadiff < 0.05) = 0;

%% downsample datadiff
downS = 5;
Ldatadiff = floor(size(datadiff,1)/downS)*downS;
datadiff_down = squeeze(mean(reshape(datadiff(1:Ldatadiff,:,:), [downS floor(size(datadiff,1)/downS) size(datadiff,2) size(datadiff,3)]),1));

%%
for i = 1:length(drivencells)
figure;
% 
%     subplot(5,5,i)
plot(tsmovavg(squeeze(datadiff_down(:,(drivencells(i)),V_cycInd)),'s',3,1),'g')
alpha(0.25)
hold on
plot(tsmovavg(squeeze(datadiff_down(:,(drivencells(i)),AV_cycInd)),'s',3,1),'k')
% hold on
vline(30,'k') 
    for iL = 1:cycles(icyc)-1
        L = (iL*(cycTime/downS))+3;
        vline(L,'k:');
        hold on
    end
hold on
plot(mean(tsmovavg(squeeze(datadiff_down(:,(drivencells(i)),V_cycInd)),'s',3,1),2),'g','LineWidth',3)
hold on
plot(mean(tsmovavg(squeeze(datadiff_down(:,(drivencells(i)),AV_cycInd)),'s',3,1),2),'k','LineWidth',3)
hold on
end

%% find signal correlation of datadiff_down (within a cell between trials)

for idrcell = 1:length(drivencells)
cell1 = drivencells(idrcell);

V_cellAvgResp = squeeze(datadiff_down(:,cell1,V_cycInd));
AV_cellAvgResp = squeeze(datadiff_down(:,cell1,AV_cycInd));

rS_V = corrcoef(V_cellAvgResp);
rS_AV = corrcoef(AV_cellAvgResp);

rS_V = tril(rS_V,-1);
rS_AV = tril(rS_AV,-1);

avgrS_V = nanmean(nanmean(rS_V(rS_V  ~= 0),2),1);
avgrS_AV = nanmean(nanmean(rS_AV(rS_AV  ~= 0),2),1);

figure;
colormap(brewermap([],'*RdBu'))
subplot(1,2,1);
imagesc(rS_V)
hold on
title(['rS_V avg corr = ' num2str(avgrS_V)])
hold on
colorbar 
caxis([-1 1]);
axis('square')
hold on
subplot(1,2,2);
imagesc(rS_AV)
title(['rS_A_V avg corr = ' num2str(avgrS_AV)])
hold on
colorbar
caxis([-1 1]);
axis('square')
end

%% calc and plot average signal correlation for a cell of interest, V vs AV
downS = 1;
cell1 = drivencells(2);
avgrS_V = zeros(length(cycles),1);
avgrS_AV = zeros(length(cycles),1);

for icyc = 1:length(cycles)
    tempdata = cycDataDFoverF_cmlvNoTarget{icyc};
    
    tempOutcome = cycTrialOutcome{icyc};
    tempVind = intersect(cycV_ind{icyc},find(strcmp(tempOutcome,'success') == 1));
    tempAVind = intersect(cycAV_ind{icyc},find(strcmp(tempOutcome,'success') == 1));
%     tempVind = cycV_ind{icyc};
%     tempAVind = cycAV_ind{icyc};    
    
    tempdatadiff = zeros(size(tempdata,1)-1,size(tempdata,3));
    tempdatadiff = diff(squeeze(tempdata(:,cell1,:)));
   
%     tempdatadiff(tempdatadiff < 0.05) = 0;
    
    tempLdatadiff = floor(size(tempdatadiff,1)/downS)*downS;
    tempdatadiff_down = squeeze(mean(reshape(tempdatadiff(1:tempLdatadiff,:,:), [downS floor(size(tempdatadiff,1)/downS) size(tempdatadiff,2) size(tempdatadiff,3)]),1));
    
    tempV_cellAvgResp = squeeze(tempdatadiff_down(:,tempVind));
    tempAV_cellAvgResp = squeeze(tempdatadiff_down(:,tempAVind));

    temprS_V = corrcoef(tempV_cellAvgResp);
    temprS_AV = corrcoef(tempAV_cellAvgResp);

    temprS_V = tril(temprS_V,-1);
    temprS_AV = tril(temprS_AV,-1);

    avgrS_V(icyc,:) = nanmean(nanmean(temprS_V(temprS_V  ~= 0),2),1);
    avgrS_AV(icyc,:) = nanmean(nanmean(temprS_AV(temprS_AV  ~= 0),2),1);
    
    steS_V(icyc,:) = (std(temprS_V(temprS_V(:)  ~= 0)))/(sqrt(sum(~isnan(temprS_V(temprS_V(:)  ~= 0)))));
    steS_AV(icyc,:) = (std(temprS_AV(temprS_AV(:)  ~= 0)))/(sqrt(sum(~isnan(temprS_AV(temprS_AV(:)  ~= 0)))));
end

figure;
errorbar(avgrS_V,steS_V,'g')
hold on
errorbar(avgrS_AV,steS_AV,'k')
%% plot this for a group of driven cells
downS = 2;
for idrcell = 1:length(drivencells) 
    
avgrS_V = zeros(length(cycles),1);
avgrS_AV = zeros(length(cycles),1);

for icyc = 1:length(cycles)
    tempdata = cycDataDFoverF_cmlvNoTarget{icyc};
    tempdata = tempdata(end-29:end,:,:);
    
    tempOutcome = cycTrialOutcome{icyc};
    tempVind = intersect(cycV_ind{icyc},find(strcmp(tempOutcome,'success') == 1));
    tempAVind = intersect(cycAV_ind{icyc},find(strcmp(tempOutcome,'success') == 1));
%     tempVind = cycV_ind{icyc};
%     tempAVind = cycAV_ind{icyc}; 
    

    tempdatadiff = zeros(size(tempdata,1)-1,size(tempdata,3));
    tempdatadiff = diff(squeeze(tempdata(:,drivencells(idrcell),:)));   
    
%     tempdatadiff(tempdatadiff < 0.05) = 0;
    
    tempLdatadiff = floor(size(tempdatadiff,1)/downS)*downS;
    tempdatadiff_down = squeeze(mean(reshape(tempdatadiff(1:tempLdatadiff,:,:), [downS floor(size(tempdatadiff,1)/downS) size(tempdatadiff,2) size(tempdatadiff,3)]),1));

    tempV_cellAvgResp = squeeze(tempdatadiff_down(:,tempVind));
    tempAV_cellAvgResp = squeeze(tempdatadiff_down(:,tempAVind));

    temprS_V = corrcoef(tempV_cellAvgResp);
    temprS_AV = corrcoef(tempAV_cellAvgResp);

    temprS_V = tril(temprS_V,-1);
    temprS_AV = tril(temprS_AV,-1);

    avgrS_V(icyc,:) = nanmean(nanmean(temprS_V(temprS_V  ~= 0),2),1);
    avgrS_AV(icyc,:) = nanmean(nanmean(temprS_AV(temprS_AV  ~= 0),2),1);
    
    steS_V(icyc,:) = (nanstd(temprS_V(temprS_V(:)  ~= 0)))/(sqrt(sum(~isnan(temprS_V(temprS_V(:)  ~= 0)))));
    steS_AV(icyc,:) = (nanstd(temprS_AV(temprS_AV(:)  ~= 0)))/(sqrt(sum(~isnan(temprS_AV(temprS_AV(:)  ~= 0)))));
end

figure;
errorbar(avgrS_V,steS_V,'g')
hold on
errorbar(avgrS_AV,steS_AV,'k')
end

%% plot this for a group of cells
downS = 1;
cellgroup = 1:nCells;
cellgroupname = 'all';
mV = zeros(length(cellgroup),1);
mAV = zeros(length(cellgroup),1);
for idrcell = 1:length(cellgroup) 
    
avgrS_V = zeros(length(cycles),1);
avgrS_AV = zeros(length(cycles),1);

for icyc = 1:length(cycles)
    tempdata = cycDataDFoverF_cmlvNoTarget{icyc};
    tempdata = tempdata(end-29:end,:,:);
    
    tempOutcome = cycTrialOutcome{icyc};
    tempVind = intersect(cycV_ind{icyc},find(strcmp(tempOutcome,'success') == 1));
    tempAVind = intersect(cycAV_ind{icyc},find(strcmp(tempOutcome,'success') == 1));
%     tempVind = cycV_ind{icyc};
%     tempAVind = cycAV_ind{icyc}; 
    

    tempdatadiff = zeros(size(tempdata,1)-1,size(tempdata,3));
    tempdatadiff = diff(squeeze(tempdata(:,cellgroup(idrcell),:)));   
    
%     tempdatadiff(tempdatadiff < 0.05) = 0;
    
    tempLdatadiff = floor(size(tempdatadiff,1)/downS)*downS;
    tempdatadiff_down = squeeze(mean(reshape(tempdatadiff(1:tempLdatadiff,:,:), [downS floor(size(tempdatadiff,1)/downS) size(tempdatadiff,2) size(tempdatadiff,3)]),1));

    tempV_cellAvgResp = squeeze(tempdatadiff_down(:,tempVind));
    tempAV_cellAvgResp = squeeze(tempdatadiff_down(:,tempAVind));

    temprS_V = corrcoef(tempV_cellAvgResp);
    temprS_AV = corrcoef(tempAV_cellAvgResp);

    temprS_V = tril(temprS_V,-1);
    temprS_AV = tril(temprS_AV,-1);

    avgrS_V(icyc,:) = nanmean(nanmean(temprS_V(temprS_V  ~= 0),2),1);
    avgrS_AV(icyc,:) = nanmean(nanmean(temprS_AV(temprS_AV  ~= 0),2),1);
    
    steS_V(icyc,:) = (nanstd(temprS_V(temprS_V(:)  ~= 0)))/(sqrt(sum(~isnan(temprS_V(temprS_V(:)  ~= 0)))));
    steS_AV(icyc,:) = (nanstd(temprS_AV(temprS_AV(:)  ~= 0)))/(sqrt(sum(~isnan(temprS_AV(temprS_AV(:)  ~= 0)))));
end
    % find slope of min cyc to max cyc correlation
mV(idrcell,:) = (avgrS_V(maxCyclesOn,:) - avgrS_V(minCyclesOn,:))/(maxCyclesOn - minCyclesOn);
mAV(idrcell,:) = (avgrS_AV(maxCyclesOn,:) - avgrS_AV(minCyclesOn,:))/(maxCyclesOn - minCyclesOn);

% figure;
% errorbar(avgrS_V,steS_V,'g')
% hold on
% errorbar(avgrS_AV,steS_AV,'k')

end

% plot difference in slope bw V and AV trials
[m_subsort m_ind] = sort(abs(mV)-abs(mAV));
figure;
scatter(1:length(cellgroup),mV(m_ind),'go')
hold on
scatter(1:length(cellgroup),mAV(m_ind),'ko')
hold on
plot(m_subsort,'r')
hold on
vline(length(cellgroup)/2)
hline(0,'k')
title([mouse ' ' date; 'cells = ' cellgroupname])
