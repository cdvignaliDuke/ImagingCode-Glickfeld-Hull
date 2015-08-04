%% find cells that respond to first stimulus, visual and auditory conditions
dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_ind = cycV_ind{1};
AV_ind = cycAV_ind{1};

preStimResp_V = zeros(size(V_ind,2),size(dataTrialStart,2));
for itrial =1:size(V_ind,1);
    for icell = 1:size(dataTrialStart,2)
        preStimResp_V(itrial,icell) = mean(dataTrialStart(1:30,icell,V_ind(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(V_ind,2),size(dataTrialStart,2));
for itrial = 1:size(V_ind,1);
    for icell = 1:size(dataTrialStart,2)
        baselineStimResp_V(itrial,icell) = mean(dataTrialStart(36:40,icell,V_ind(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);

preStimResp_AV = zeros(size(AV_ind,2),size(dataTrialStart,2));
for itrial =1:size(AV_ind,1);
    for icell = 1:size(dataTrialStart,2)
        preStimResp_AV(itrial,icell) = mean(dataTrialStart(1:30,icell,AV_ind(itrial)),1);
    end
end

baselineStimResp_AV = zeros(size(AV_ind,2),size(dataTrialStart,2));
for itrial = 1:size(AV_ind,1);
    for icell = 1:size(dataTrialStart,2)
        baselineStimResp_AV(itrial,icell) = mean(dataTrialStart(36:40,icell,AV_ind(itrial)),1);
    end
end

baselineStimRespTtest_AV= ttest(preStimResp_AV,baselineStimResp_AV,'alpha', 0.01);
baselineStimRespIndex_AV = find(baselineStimRespTtest_AV == 1);

%% tuning
dirFolder = '005';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
cd(fileSave);
% load('oriTuningPreferences.mat')
load('TuningPreferences.mat')

cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsSelectZero_dir = intersect(dirSlctvCells,cellsPrefZero);
cellsSelectZero_ori = intersect(oriSlctvCells,cellsPrefZero);
cellsSelectZero = union(cellsSelectZero_dir,cellsSelectZero_ori);
% respCellsSelectZero = intersect(baselineStimRespIndex_V,cellsSelectZero);

cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
cellsSelectNinety_dir = intersect(dirSlctvCells,cellsPrefNinety);
cellsSelectNinety_ori = intersect(oriSlctvCells,cellsPrefNinety);
cellsSelectNinety = union(cellsSelectNinety_dir,cellsSelectNinety_ori);

%%
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,intersect(cellsPrefNinety,baselineStimRespIndex_V),V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,:,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,intersect(cellsPrefNinety,baselineStimRespIndex_V),AV_cycInd),3),2);
    errbar_V = (std(mean(dataDFoverF(:,intersect(cellsPrefNinety,baselineStimRespIndex_V),V_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,intersect(cellsPrefNinety,baselineStimRespIndex_V),V_cycInd),3)));
%     errbar_A = (std(mean(dataDFoverF(:,:,A_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,:,A_cycInd),3)));
    errbar_AV = (std(mean(dataDFoverF(:,intersect(cellsPrefNinety,baselineStimRespIndex_V),AV_cycInd),2),[],3))/(sqrt(size(dataDFoverF(:,intersect(cellsPrefNinety,baselineStimRespIndex_V),AV_cycInd),3)));
    subplot(3,3,icyc);
%     plot(V_avg(20:end,:),'g');
%     hold on
% %     plot(A_avg(20:end,:),'r');
%     hold on
%     plot(AV_avg(20:end,:),'m');
%     hold on
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g')
    hold on
%     errorbar(A_avg(20:end,:),errbar_A(20:end,:),'r')
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
        title([num2str(size(dataDFoverF(:,intersect(baselineStimRespIndex_V,cellsSelectZero),V_cycInd),2)) ' cells'])
    else
%     title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(A_cycInd)) ' auditory trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
    title([num2str(length(V_cycInd)) ' visual trials; ' num2str(length(AV_cycInd)) ' vis+aud trials'])
    end
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
end