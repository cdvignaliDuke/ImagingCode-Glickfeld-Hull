edit FlashingStim_dataSortedByCycle_combineDatasetsSameType.m
edit successTrials.m

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('cellSelectivityIndices.mat')

cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 7);
cellsPrefNinety = find(dirPref_ind == 4 | dirPref_ind == 10);
cellsSlctvZero = intersect(cellsPrefZero,dirSlctvCells);


%% plot cell traces

dataCellsTrials4cycles = cycDataDFoverF{4};


cellsperpage = 11;
pages = ceil(size(dataCellsTrials4cycles,2)/cellsperpage);
cycles = 4;

cell = 1;
for ifig = 1:pages
figure;
for iplot = 1:cellsperpage
    subplot(3,4,iplot)
    trials = find(tCyclesOn == cycles);
    V_cycInd = find(ismember(trials,V_ind));
    plot(mean(dataCellsTrials4cycles(:,cell,V_cycInd),3),'g')
    hold on
    AV_cycInd = find(ismember(trials,AV_ind));
    plot(mean(dataCellsTrials4cycles(:,cell,AV_cycInd),3),'m')
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles*cycTime+30),'c');
    xlim([0 size(dataCellsTrials4cycles,1)])
    title(['Cell ' num2str(cell)])
    if iplot == 1
        title([num2str(length(V_cycInd)) ' vis, ' num2str(length(AV_cycInd)) 'aud+vis'])
    end
    if ismember(cell,cellsPrefZero) > 0
        set(subplot(3,4,iplot),'color',[0.9 0.9 0.9])
    end
    cell = cell+1;
end    
end

dataCellsTrials4cycles_notarget = cycDataDFoverF_cmlvNoTarget{4};
V4_ind = cycV_ind{4};
AV4_ind = cycAV_ind{4};

cell = 1;
for ifig = 1:pages
figure;
for iplot = 1:cellsperpage
    subplot(3,4,iplot)
    plot(mean(dataCellsTrials4cycles_notarget(:,cell,V4_ind),3),'g')
    hold on
    plot(mean(dataCellsTrials4cycles_notarget(:,cell,AV4_ind),3),'m')
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles*cycTime+30),'c');
    xlim([0 size(dataCellsTrials4cycles_notarget,1)])
    title(['Cell ' num2str(cell)])
    if iplot == 1
        title([num2str(length(V_ind)) ' vis, ' num2str(length(AV_ind)) 'aud+vis'])
    end
    if ismember(cell,cellsPrefZero) > 0
        set(subplot(3,4,iplot),'color',[0.9 0.9 0.9])
    end
    cell = cell+1;
end    
end
