edit('FlashingStim_dataSortedByCycle_combineDatasets.m')

DirFolder = '005';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, DirFolder);
cd(fileSave);
% load('cellSelectivityIndices.mat')
load('TuningPreferences.mat')

cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);

cellsSelectZero_dir = intersect(dirSlctvCells,cellsPrefZero);
cellsSelectZero_ori = intersect(oriSlctvCells,cellsPrefZero);
cellsSelectZero = union(cellsSelectZero_dir,cellsSelectZero_ori);

cellsSelectNinety_dir = intersect(dirSlctvCells,cellsPrefNinety);
cellsSelectNinety_ori = intersect(oriSlctvCells,cellsPrefNinety);
cellsSelectNinety = union(cellsSelectNinety_dir,cellsSelectNinety_ori);

% respCellsSelectZero = intersect(baselineStimRespIndex_V,cellsSelectZero);

%% plot dF/F for specific sets of cells (cumulative, no target)
% one specific cell
cell = 116;
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(dataDFoverF(:,cell,V_cycInd),3);
    A_avg = mean(dataDFoverF(:,cell,A_cycInd),3);
    AV_avg = mean(dataDFoverF(:,cell,AV_cycInd),3);
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_cycInd)) ' vis trials; ' num2str(length(A_cycInd)) ' aud trials ' num2str(length(AV_cycInd)) ' vis+aud trials']);
    hold on
end

% plot cells that are orientation or direction selective, prefer zero
% degrees orientation, and significantly respond to baseline stim
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
%     V_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,A_cycInd),3),2);
%     AV_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,AV_cycInd),3),2);
    V_avg = mean(mean(dataDFoverF(:,cellsSelectZero,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,cellsSelectZero,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cellsSelectZero,AV_cycInd),3),2);
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(AV_ind)) ' auditory trials'])
    hold on
    if icyc == 1
        title(['n = ' num2str(length(cellsSelectZero)) '; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(AV_ind)) ' auditory trials']);
    end
end

figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
%     A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
%     V_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,A_cycInd),3),2);
%     AV_avg = mean(mean(dataDFoverF(:,respCellsSelectZero,AV_cycInd),3),2);
    V_avg = mean(mean(dataDFoverF(:,cellsSelectNinety,V_cycInd),3),2);
%     A_avg = mean(mean(dataDFoverF(:,cellsSelectZero,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,cellsSelectNinety,AV_cycInd),3),2);
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(AV_ind)) ' auditory trials'])
    hold on
    if icyc == 1
        title(['n = ' num2str(length(cellsSelectNinety)) '; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(AV_ind)) ' auditory trials']);
    end
end
%% response reliability
figure;
dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_cycInd = cycV_ind{1};
% A_cycInd = cycA_ind{1};
AV_cycInd = cycAV_ind{1};
% V_avg = mean(mean(dataTrialStart(:,respCellsSelectZero,V_cycInd),3),2);
% A_avg = mean(mean(dataTrialStart(:,respCellsSelectZero,A_cycInd),3),2);
% AV_avg = mean(mean(dataTrialStart(:,respCellsSelectZero,AV_cycInd),3),2);
V_avg = mean(mean(dataTrialStart(:,cellsSelectZero,V_cycInd),3),2);
% A_avg = mean(mean(dataTrialStart(:,cellsSelectZero,A_cycInd),3),2);
AV_avg = mean(mean(dataTrialStart(:,cellsSelectZero,AV_cycInd),3),2);
plot(V_avg,'g');
hold on
% plot(A_avg,'r');
hold on
plot(AV_avg,'m');
hold on
vline(30,'k')
hold on
for i = 1:cycles(1)-1
    L = (i*cycTime)+30;
    vline(L,'k:');
    hold on
end
vline((cycles(1)*cycTime+30),'c');
hold on
% title([num2str(cycles(1)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials'])
hold on

%plot cell averages
% plotRows = ceil(size(respCellsSelectZero,1)/4);
% 
% V_cellAvg = zeros(size(dataTrialStart,1),length(respCellsSelectZero)); 
% A_cellAvg = zeros(size(dataTrialStart,1),length(respCellsSelectZero));
% AV_cellAvg = zeros(size(dataTrialStart,1),length(respCellsSelectZero));
% for i = 1:size(respCellsSelectZero,1)
%     cell = respCellsSelectZero(i);
%     V_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,V_cycInd),3));
%     A_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,A_cycInd),3));
%     AV_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,AV_cycInd),3));
% end
% 
% figure;
% for i = 1:size(respCellsSelectZero,1)
%     cell = respCellsSelectZero(i)
%     subplot(4,plotRows,i)
%     plot(squeeze(V_cellAvg(:,i,:)),'g');
%     hold on
%     plot(squeeze(A_cellAvg(:,i,:)),'r');
%     hold on
%     plot(squeeze(AV_cellAvg(:,i,:)),'m');
%     hold on
%     vline(30,'k')
%     hold on
%     for i = 1:cycles(1)-1
%         L = (i*cycTime)+30;
%         vline(L,'k:');
%         hold on
%     end
%     vline((cycles(1)*cycTime+30),'c');
%     hold on
%     title(['Cell ' num2str(cell)]);
%     hold on
% end


%% cells prefer 0 deg
plotRows = ceil(size(cellsSelectZero,1)/4);
dataTrialStart = cycDataDFoverF_cmlvNoTarget{7};
V_cycInd = cycV_ind{7};
% A_cycInd = cycA_ind{end};
AV_cycInd = cycAV_ind{7};

V_cellAvg = zeros(size(dataTrialStart,1),length(cellsSelectZero)); 
% A_cellAvg = zeros(size(dataTrialStart,1),length(cellsSelectZero));
AV_cellAvg = zeros(size(dataTrialStart,1),length(cellsSelectZero));

for i = 1:size(cellsSelectZero)
    cell = cellsSelectZero(i);
    V_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,V_cycInd),3));
    AV_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,AV_cycInd),3));
end


figure;
for i = 1:size(cellsSelectZero,1)
    cell = cellsSelectZero(i);
    subplot(4,plotRows,i)
    plot(squeeze(V_cellAvg(:,i,:)),'g');
    hold on
%     plot(squeeze(A_cellAvg(:,i,:)),'r');
    hold on
    plot(squeeze(AV_cellAvg(:,i,:)),'m');
    hold on
    vline(30,'k')
    hold on
    for i = 1:cycles(7)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(7)*cycTime+30),'c');
    hold on
    title(['Cell ' num2str(cell)]);
    hold on
end


%calculate the standard devation over timecourse of trial start for each
%resposive, orientation selective, zero prefering cell

V_std = std(dataTrialStart(:,:,V_cycInd),0,3);
V_cycInd = cycV_ind{7};
% A_cycInd = cycA_ind{7};
AV_cycInd = cycAV_ind{7};
V_std = std(dataTrialStart(:,cellsSelectZero,V_cycInd),0,3);
% A_std = std(dataTrialStart(:,cellsSelectZero,A_cycInd),0,3);
AV_std = std(dataTrialStart(:,cellsSelectZero,AV_cycInd),0,3);

figure;
for i = 1:size(cellsSelectZero,1)
    cell = cellsSelectZero(i)
    subplot(4,plotRows,i)
    plot(squeeze(V_std(:,i,:)),'g');
    hold on
%     plot(squeeze(A_std(:,i,:)),'r');
    hold on
    plot(squeeze(AV_std(:,i,:)),'m');
    hold on
    vline(30,'k')
    hold on
    for i = 1:cycles(7)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(7)*cycTime+30),'c');
    hold on
    title(['Cell ' num2str(cell)]);
    hold on
end

%calculate the coeffient of variation over timecourse of trial start for each
%resposive, orientation selective, zero prefering cell

V_CV = zeros(size(V_std));
% A_CV = zeros(size(A_std));
AV_CV = zeros(size(AV_std));
for i = 1:size(cellsSelectZero,1)    
    V_CV(:,i) = bsxfun(@rdivide,V_std(:,i),V_cellAvg(:,i));
%     A_CV(:,i) = bsxfun(@rdivide,A_std(:,i),A_cellAvg(:,i));
    AV_CV(:,i) = bsxfun(@rdivide,AV_std(:,i),AV_cellAvg(:,i));
end

figure;
for i = 1:size(cellsSelectZero,1)
    cell = cellsSelectZero(i)
    subplot(4,plotRows,i)
    plot(squeeze(V_CV(:,i,:)),'g');
    hold on
%     plot(squeeze(A_CV(:,i,:)),'r');
    hold on
    plot(squeeze(AV_CV(:,i,:)),'m');
    hold on
    vline(30,'k')
    hold on
    for i = 1:cycles(7)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(7)*cycTime+30),'c');
    hold on
    title(['Cell ' num2str(cell)]);
    hold on
end

%% cells prefer 90 deg
%% cells prefer 0 deg
plotRows = ceil(size(cellsSelectNinety(21:40),1)/4);
dataTrialStart = cycDataDFoverF_cmlvNoTarget{7};
V_cycInd = cycV_ind{7};
% A_cycInd = cycA_ind{end};
AV_cycInd = cycAV_ind{7};

V_cellAvg = zeros(size(dataTrialStart,1),length(cellsSelectNinety)); 
% A_cellAvg = zeros(size(dataTrialStart,1),length(cellsSelectNinety));
AV_cellAvg = zeros(size(dataTrialStart,1),length(cellsSelectNinety));

for i = 1:size(cellsSelectNinety)
    cell = cellsSelectNinety(i);
    V_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,V_cycInd),3));
    AV_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,AV_cycInd),3));
end


figure;
for i = 1:size(cellsSelectNinety(21:40),1)
    cell = cellsSelectNinety(i);
    subplot(4,plotRows,i)
    plot(squeeze(V_cellAvg(:,i,:)),'g');
    hold on
%     plot(squeeze(A_cellAvg(:,i,:)),'r');
    hold on
    plot(squeeze(AV_cellAvg(:,i,:)),'m');
    hold on
    vline(30,'k')
    hold on
    for i = 1:cycles(7)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(7)*cycTime+30),'c');
    hold on
    title(['Cell ' num2str(cell)]);
    hold on
end


%calculate the standard devation over timecourse of trial start for each
%resposive, orientation selective, zero prefering cell

V_std = std(dataTrialStart(:,:,V_cycInd),0,3);
V_cycInd = cycV_ind{7};
% A_cycInd = cycA_ind{7};
AV_cycInd = cycAV_ind{7};
V_std = std(dataTrialStart(:,cellsSelectNinety,V_cycInd),0,3);
% A_std = std(dataTrialStart(:,cellsSelectZero,A_cycInd),0,3);
AV_std = std(dataTrialStart(:,cellsSelectNinety,AV_cycInd),0,3);

figure;
for i = 1:size(cellsSelectNinety,1)
    cell = cellsSelectNinety(i);
    subplot(4,plotRows,i)
    plot(squeeze(V_std(:,i,:)),'g');
    hold on
%     plot(squeeze(A_std(:,i,:)),'r');
    hold on
    plot(squeeze(AV_std(:,i,:)),'m');
    hold on
    vline(30,'k')
    hold on
    for i = 1:cycles(7)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(7)*cycTime+30),'c');
    hold on
    title(['Cell ' num2str(cell)]);
    hold on
end