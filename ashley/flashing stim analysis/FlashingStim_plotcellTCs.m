edit FlashingStim_dataSortedByCycle_combineDatasetsSameType.m

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('cellSelectivityIndices.mat')
% 
% cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
% cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
% cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 7);
% cellsPrefNinety = find(dirPref_ind == 4 | dirPref_ind == 10);
cellsSlctvZero = intersect(cellsPrefZero,dirSlctvCells);


%% plot cell traces
cycTrialLength = 6;
cyc = find(cycles == cycTrialLength);
dataCellsTrials4cycles = cycDataDFoverF{cyc};


cellsperpage = 15;
pages = ceil(size(dataCellsTrials4cycles,2)/cellsperpage);


cell = 1;
for ifig = 1:pages
figure;
for iplot = 1:32
    subplot(8,4,iplot)
    trials = find(tCyclesOn == cyc);
    V_cycInd = find(ismember(trials,V_ind));
    plot(mean(dataCellsTrials4cycles(:,cell,V_cycInd),3),'g')
    hold on
%     plot(squeeze(dataCellsTrials4cycles(:,cell,V_cycInd)),'g');
%     hold on
    AV_cycInd = find(ismember(trials,AV_ind));
    plot(mean(dataCellsTrials4cycles(:,cell,AV_cycInd),3),'m')
    hold on
%     plot(squeeze(dataCellsTrials4cycles(:,cell,AV_cycInd)),'m');
%     hold on
    vline(30,'k');
    hold on
    for i = 1:cyc-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cyc*cycTime+30),'c');
    xlim([0 size(dataCellsTrials4cycles,1)])
    title(['Cell ' num2str(cell)])
    if iplot == 1
        title([num2str(length(V_cycInd)) ' vis, ' num2str(length(AV_cycInd)) 'aud+vis'])
    end
    if ismember(cell,cellsPrefZero) > 0
        set(subplot(8,4,iplot),'color',[0.9 0.9 0.9])
    end
    if ismember(cell,cellsPrefNinety) > 0
        set(subplot(8,4,iplot),'color',[0.8 0.8 0.8])
    end
    cell = cell+1;
end    
end

dataCellsTrials4cycles_notarget = cycDataDFoverF_cmlvNoTarget{cyc};
V4_ind = cycV_ind{cyc};
AV4_ind = cycAV_ind{cyc};

cell = 1;
for ifig = 1:pages
figure;
for iplot = 1:cellsperpage
    subplot(4,4,iplot)
    plot(mean(dataCellsTrials4cycles_notarget(:,cell,V4_ind),3),'g')
    hold on
    plot(mean(dataCellsTrials4cycles_notarget(:,cell,AV4_ind),3),'m')
    hold on
    vline(30,'k');
    hold on
    for i = 1:cyc-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline(((cyc)*cycTime+30),'c');
    xlim([0 size(dataCellsTrials4cycles_notarget,1)+5])
    ylim([-0.05 0.1])
    title(['Cell ' num2str(cell)])
    if iplot == 1
        title([num2str(length(V_ind)) ' vis, ' num2str(length(AV_ind)) 'aud+vis'])
    end
%     if ismember(cell,cellsPrefZero) > 0
%         set(subplot(4,4,iplot),'color',[0.9 0.9 0.9])
%     end
%     if ismember(cell,cellsPrefNinety) > 0
%         set(subplot(4,4,iplot),'color',[0.8 0.8 0.8])
%     end
    cell = cell+1;
end    
end

%% select cells
cells = [ 9 18 29 36 40 41 44 50 55];
cellstart = 1;
for ifig = 1:pages
figure;
for iplot = 1:9
    cell = cells(cellstart);
    subplot(3,3,iplot)
    plot(mean(dataCellsTrials4cycles_notarget(:,cell,V4_ind),3),'g')
    hold on
    plot(mean(dataCellsTrials4cycles_notarget(:,cell,AV4_ind),3),'m')
    hold on
    vline(30,'k');
    hold on
    for i = 1:cyc
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline(((cyc+1)*cycTime+30),'c');
    xlim([0 size(dataCellsTrials4cycles_notarget,1)+5])
    ylim([-0.05 0.1])
    title(['Cell ' num2str(cell)])
    if iplot == 1
        title([num2str(length(V_ind)) ' vis, ' num2str(length(AV_ind)) 'aud+vis; ' num2str(cycTrialLength) ' cycles'])
    end
    if ismember(cell,cellsPrefZero) > 0
        set(subplot(3,3,iplot),'color',[0.9 0.9 0.9])
    end
    if ismember(cell,cellsPrefNinety) > 0
        set(subplot(3,3,iplot),'color',[0.8 0.8 0.8])
    end
    cellstart = cellstart+1;
end    
end
