edit('FlashingStim_dataSortedByCycle_combineDatasets.m')
% edit FlashingStim_dataSortedByCycle_combineDatasetsSameType.m

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);

dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_ind = cycV_ind{1};
AV_ind = cycAV_ind{1};
%% find cells that respond to first stimulus, visual and auditory conditions
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

baselineStimRespTtest_AV = ttest(preStimResp_AV,baselineStimResp_AV,'alpha', 0.01);
baselineStimRespIndex_AV = find(baselineStimRespTtest_AV == 1);

%% find cells with 0 deg orientation preference
% dirTuningFolder = '005';
% CD = ['Z:\analysis\' mouse '\two-photon imaging\' date '\' dirTuningFolder];
% cd(CD);
% load('oriTuningPreferences.mat')

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);

load('cellSelectivityIndices.mat')

cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
cellsSlctvZero = intersect(cellsPrefZero,dirSlctvCells);
% cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 7);
% cellsPrefNinety = find(dirPref_ind == 4 | dirPref_ind == 10);
% cellsSlctvZero = intersect(cellsPrefZero,dirSlctvCells);

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% save('cellSelectivityIndices.mat', 'baselineStimRespIndex_V','baselineStimRespIndex_AV','dirSlctvCells','dirPref_ind','oriSlctvCells','oriPref_ind')

%% plot cycle-type traces for specific cells
cell = 10;
figure;
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = find(tCyclesOn == cycles(icyc));
    V_cycInd = find(ismember(trials,V_ind));
%     A_cycInd = find(ismember(trials,A_ind));
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(DataDFoverF(:,cell,V_cycInd),3);
%     A_avg = mean(DataDFoverF(:,cell,A_cycInd),3);
    AV_avg = mean(DataDFoverF(:,cell,AV_cycInd),3);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
%     plot(A_avg,'r');
    hold on
    plot(AV_avg,'m');
    hold on
%     title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
end

cell = 10;
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(dataDFoverF(:,cell,V_cycInd),3);
%     A_avg = mean(dataDFoverF(:,cell,A_cycInd),3);
    AV_avg = mean(dataDFoverF(:,cell,AV_cycInd),3);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
%     plot(A_avg,'r');
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
%     title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
end

%% plot averages of visually responsive cells
figure;
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = find(tCyclesOn == cycles(icyc));
    V_cycInd = find(ismember(trials,V_ind));
    A_cycInd = find(ismember(trials,A_ind));
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(mean(DataDFoverF(:,baselineStimRespIndex_V,V_cycInd),3),2);
    A_avg = mean(mean(DataDFoverF(:,baselineStimRespIndex_V,A_cycInd),3),2);
    AV_avg = mean(mean(DataDFoverF(:,baselineStimRespIndex_V,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(A_avg,'r');
    hold on
    plot(AV_avg,'m');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,baselineStimRespIndex_AV,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,baselineStimRespIndex_V,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,baselineStimRespIndex_AV,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
%     plot(A_avg,'r');
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
%     title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,baselineStimRespIndex_A,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,baselineStimRespIndex_V,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,baselineStimRespIndex_A,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
%     plot(A_avg,'r');
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
%     title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
end
%% plot averages of cell with 0 deg ori pref
figure;
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = find(tCyclesOn == cycles(icyc));
    V_cycInd = find(ismember(trials,V_ind));
    A_cycInd = find(ismember(trials,A_ind));
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(mean(DataDFoverF(:,cellsPrefZero,V_cycInd),3),2);
    A_avg = mean(mean(DataDFoverF(:,cellsPrefZero,A_cycInd),3),2);
    AV_avg = mean(mean(DataDFoverF(:,cellsPrefZero,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(A_avg,'r');
    hold on
    plot(AV_avg,'m');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,cellsPrefZero,V_cycInd),3),2);
    A_avg = mean(mean(dataDFoverF(:,cellsPrefZero,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cellsPrefZero,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg(20:end,:),'g');
    hold on
    plot(A_avg(20:end,:),'r');
    hold on
    plot(AV_avg(20:end,:),'m');
    hold on
    vline(10,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+10;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+10),'c');
    hold on
    title(['n = ' num2str(length(cellsPrefZero)) '; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
end

%% plot averages of cell with 90 deg ori pref
figure;
for icyc = 1:length(cycles)
    DataDFoverF = cycDataDFoverF{icyc};
    trials = find(tCyclesOn == cycles(icyc));
    V_cycInd = find(ismember(trials,V_ind));
    A_cycInd = find(ismember(trials,A_ind));
    AV_cycInd = find(ismember(trials,AV_ind));
    V_avg = mean(mean(DataDFoverF(:,cellsPrefNinety,V_cycInd),3),2);
    A_avg = mean(mean(DataDFoverF(:,cellsPrefNinety,A_cycInd),3),2);
    AV_avg = mean(mean(DataDFoverF(:,cellsPrefNinety,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg,'g');
    hold on
    plot(A_avg,'r');
    hold on
    plot(AV_avg,'m');
    hold on
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,cellsPrefNinety,V_cycInd),3),2);
    A_avg = mean(mean(dataDFoverF(:,cellsPrefNinety,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cellsPrefNinety,AV_cycInd),3),2);
    subplot(3,3,icyc);
    plot(V_avg(20:end,:),'g');
    hold on
    plot(A_avg(20:end,:),'r');
    hold on
    plot(AV_avg(20:end,:),'m');
    hold on
    vline(10,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+10;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+10),'c');
    hold on
    title(['n = ' num2str(length(cellsPrefNinety)) '; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
end

%% use error bars
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    errbar_V = std(mean(dataDFoverF(:,cellsPrefZero,V_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,cellsPrefZero,V_cycInd),3));
    errbar_A = std(mean(dataDFoverF(:,cellsPrefZero,A_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,cellsPrefZero,A_cycInd),3));
    errbar_AV = std(mean(dataDFoverF(:,cellsPrefZero,AV_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,cellsPrefZero,AV_cycInd),3));
    V_avg = mean(mean(dataDFoverF(:,cellsPrefZero,V_cycInd),3),2);
    A_avg = mean(mean(dataDFoverF(:,cellsPrefZero,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cellsPrefZero,AV_cycInd),3),2);
    subplot(3,3,icyc);
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g');
    hold on
    errorbar(A_avg(20:end,:),errbar_A(20:end,:),'r');
    hold on
    errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'m');
    hold on
    vline(10,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+10;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+10),'c');
    hold on
    title(['n = ' num2str(length(cellsPrefZero)) '; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
%     title(['n = ' num2str(length(cellsPrefZero)) '; ' num2str(length(V_cycInd)) ' vis trials; '  num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    errbar_V = std(mean(dataDFoverF(:,cellsPrefNinety,V_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,cellsPrefNinety,V_cycInd),3));
    errbar_A = std(mean(dataDFoverF(:,cellsPrefNinety,A_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,cellsPrefNinety,A_cycInd),3));
    errbar_AV = std(mean(dataDFoverF(:,cellsPrefNinety,AV_cycInd),2),[],3)/sqrt(size(dataDFoverF(:,cellsPrefNinety,AV_cycInd),3));
    V_avg = mean(mean(dataDFoverF(:,cellsPrefNinety,V_cycInd),3),2);
    A_avg = mean(mean(dataDFoverF(:,cellsPrefNinety,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cellsPrefNinety,AV_cycInd),3),2);
    subplot(3,3,icyc);
    errorbar(V_avg(20:end,:),errbar_V(20:end,:),'g');
    hold on
    errorbar(A_avg(20:end,:),errbar_A(20:end,:),'r');
    hold on
    errorbar(AV_avg(20:end,:),errbar_AV(20:end,:),'m');
    vline(10,'k')
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+10;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+10),'c');
    hold on
    title(['n = ' num2str(length(cellsPrefNinety)) '; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
%     title(['n = ' num2str(length(cellsPrefNinety)) '; ' num2str(length(V_cycInd)) ' vis trials; '  num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
    xlim([0 length(V_avg(20:end,:))+5])
    ylim([-0.05 0.05])
end
