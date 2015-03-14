edit('FlashingStim_dataSortedByCycle_combineDatasets.m')

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('cellSelectivityIndices.mat')

cellsPrefNinety = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsSelectNinety_dir = intersect(dirSlctvCells,cellsPrefNinety);
cellsSelectNinety_ori = intersect(oriSlctvCells,cellsPrefNinety);
cellsSelectNinety = union(cellsSelectNinety_dir,cellsSelectNinety_ori);
respCellsSelectNinety = intersect(baselineStimRespIndex_V,cellsSelectNinety);

%% plot dF/F for specific sets of cells (cumulative, no target)
% one specific cell
cell = 40;
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

% plot cells that are orientation or direction selective, prefer ninety
% degrees orientation, and significantly respond to baseline stim
figure;
for icyc = 1:length(cycles)
    dataDFoverF = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = cycV_ind{icyc};
    A_cycInd = cycA_ind{icyc};
    AV_cycInd = cycAV_ind{icyc};
    V_avg = mean(mean(dataDFoverF(:,respCellsSelectNinety,V_cycInd),3),2);
    A_avg = mean(mean(dataDFoverF(:,respCellsSelectNinety,A_cycInd),3),2);
    AV_avg = mean(mean(dataDFoverF(:,respCellsSelectNinety,AV_cycInd),3),2);
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials'])
    hold on
end

%% plot cell averages
dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_cycInd = cycV_ind{1};
A_cycInd = cycA_ind{1};
AV_cycInd = cycAV_ind{1};

plotRows = ceil(size(respCellsSelectZero,1)/4);

V_cellAvg = zeros(size(dataTrialStart,1),length(respCellsSelectNinety)); 
A_cellAvg = zeros(size(dataTrialStart,1),length(respCellsSelectNinety));
AV_cellAvg = zeros(size(dataTrialStart,1),length(respCellsSelectNinety));
for i = 1:size(respCellsSelectNinety,1)
    cell = respCellsSelectNinety(i);
    V_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,V_cycInd),3));
    A_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,A_cycInd),3));
    AV_cellAvg(:,i) = squeeze(mean(dataTrialStart(:,cell,AV_cycInd),3));
end

figure;
for i = 1:size(respCellsSelectNinety,1)
    cell = respCellsSelectNinety(i)
    subplot(4,plotRows,i)
    plot(squeeze(V_cellAvg(:,i,:)),'g');
    hold on
    plot(squeeze(A_cellAvg(:,i,:)),'r');
    hold on
    plot(squeeze(AV_cellAvg(:,i,:)),'m');
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
    title(['Cell ' num2str(cell)]);
    hold on
end

