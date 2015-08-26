% run(FlashingStim_dataSortedByCycle_combineDatasets.m)
%%
% cells selective for 90 deg, 0 deg, driven by baseline stim, "driven
% cells", non-selective cells
%%
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' date '\PreTargetAvgTraces'];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' date]);
        mkdir('PreTargetAvgTraces')
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(date,'PreTargetAvgTraces')
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
cellgroupname = 'ori or dir selective';
figName = 'avgDFoverFpretarget_orislctv_success';
%%
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
%     V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'ignore')));
%     AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'ignore')));
%     V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'failure')));
%     AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'failure')));
    V_avg = mean(mean(dataDFoverF(:,cells,V_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cells,AV_cycInd),3),2);
    errbar_V = (std(mean(dataDFoverF(:,cells,V_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cells,V_cycInd),3)));
    errbar_AV = (std(mean(dataDFoverF(:,cells,AV_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cells,AV_cycInd),3)));
    
    subplot(3,4,icyc);
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g')
    hold on
    hold on
    errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'k')
    hold on
    vline(10,'k')
    hold on
    
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+11;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+11),'c');
    hold on
    
    if icyc == 1
        title([num2str(length(cells)) ' cells'])
    else
    title({[num2str(length(V_cycInd)) ' visual trials; ']; [num2str(length(AV_cycInd)) ' vis+aud trials']})
    end
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
%     ylim([-0.05 0.05])
end
suptitle([mouse '; ' date '; ' cellgroupname])

print([fnout ['\' figName '.pdf']], '-dpdf')
%%
% cells = drivencells;
% cellgroupname = 'driven cells';
% 
% figure;
% for icyc = 1:length(cycles)
%     dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
%     V1_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
%     V2_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'ignore')));
%     V1_avg = mean(mean(dataDFoverF(:,cells,V1_cycInd),3),2);
%     V2_avg = mean(mean(dataDFoverF(:,cells,V2_cycInd),3),2);
%     errbar_V1 = (std(mean(dataDFoverF(:,cells,V1_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cells,V1_cycInd),3)));
%     errbar_V2 = (std(mean(dataDFoverF(:,cells,V2_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cells,V2_cycInd),3)));
%     
%     subplot(3,4,icyc);
%     errorbar(V1_avg(20:end,:),errbar_V1(20:end,:),'g')
%     hold on
%     hold on
%     errorbar(V2_avg(20:end,:),errbar_V2(20:end,:),'k')
%     hold on
%     vline(10,'k')
%     hold on
%     
%     for i = 1:cycles(icyc)-1
%         L = (i*cycTime)+11;
%         vline(L,'k:');
%         hold on
%     end
%     vline((cycles(icyc)*cycTime+11),'c');
%     hold on
%     
%     if icyc == 1
%         title([num2str(size(dataDFoverF,2)) ' cells'])
%     else
%     title([num2str(length(V1_cycInd)) ' vis success trials; ' num2str(length(V2_cycInd)) ' vis ignore trials'])
%     end
%     hold on
%     xlim([0 length(V1_avg(20:end,:))+5])
% %     ylim([-0.05 0.05])
% end
% suptitle([mouse '; ' date '; ' cellgroupname])
% %%
% cells = drivencells;
% cellgroupname = 'driven cells';
% 
% figure;
% for icyc = 1:length(cycles)
%     dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
%     AV1_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
%     AV2_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'ignore')));
%     AV1_avg = mean(mean(dataDFoverF(:,cells,AV1_cycInd),3),2);
%     AV2_avg = mean(mean(dataDFoverF(:,cells,AV2_cycInd),3),2);
%     errbar_AV1 = (std(mean(dataDFoverF(:,cells,AV1_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cells,AV1_cycInd),3)));
%     errbar_AV2 = (std(mean(dataDFoverF(:,cells,AV2_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,cells,AV2_cycInd),3)));
%     
%     subplot(3,4,icyc);
%     errorbar(AV1_avg(20:end,:),errbar_AV1(20:end,:),'g')
%     hold on
%     hold on
%     errorbar(AV2_avg(20:end,:),errbar_AV2(20:end,:),'k')
%     hold on
%     vline(10,'k')
%     hold on
%     
%     for i = 1:cycles(icyc)-1
%         L = (i*cycTime)+11;
%         vline(L,'k:');
%         hold on
%     end
%     vline((cycles(icyc)*cycTime+11),'c');
%     hold on
%     
%     if icyc == 1
%         title([num2str(size(dataDFoverF,2)) ' cells'])
%     else
%     title([num2str(length(AV1_cycInd)) ' vis+aud success trials; ' num2str(length(AV2_cycInd)) ' vis+aud ignore trials'])
%     end
%     hold on
%     xlim([0 length(AV1_avg(20:end,:))+5])
% %     ylim([-0.05 0.05])
% end
% suptitle([mouse '; ' date '; ' cellgroupname])