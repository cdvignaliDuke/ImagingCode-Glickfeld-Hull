edit('FlashingStim_dataSortedByCycle_combineDatasets.m')

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('cellSelectivityIndices.mat')

cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsSelectZero_dir = intersect(dirSlctvCells,cellsPrefZero);
cellsSelectZero_ori = intersect(oriSlctvCells,cellsPrefZero);
cellsSelectZero = union(cellsSelectZero_dir,cellsSelectZero_ori);
respCellsSelectZero = intersect(baselineStimRespIndex_V,cellsSelectZero);

%% plot dF/F for specific sets of cells (cumulative, no target)
% one specific cell
cell = 121;
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
    title([num2str(cycles(icyc)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials'])
    hold on
end

%% response reliability
figure;
dataTrialStart = cycDataDFoverF_cmlvNoTarget{1};
V_cycInd = cycV_ind{1};
A_cycInd = cycA_ind{1};
AV_cycInd = cycAV_ind{1};
V_avg = mean(mean(dataTrialStart(:,respCellsSelectZero,V_cycInd),3),2);
A_avg = mean(mean(dataTrialStart(:,respCellsSelectZero,A_cycInd),3),2);
AV_avg = mean(mean(dataTrialStart(:,respCellsSelectZero,AV_cycInd),3),2);
plot(V_avg,'g');
hold on
plot(A_avg,'r');
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
title([num2str(cycles(1)) ' cycles; ' num2str(length(V_ind)) ' visual trials; ' num2str(length(A_ind)) ' auditory trials'])
hold on

%calculate the standard devation over timecourse of trial start for each
%resposive, orientation selective, zero prefering cell
dataTrialStart_std = std(dataTrialStart,0,3);



