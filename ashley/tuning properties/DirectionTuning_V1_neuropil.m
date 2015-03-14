%% get dF/F of cells from direction tuning
edit DirectionTuning_V1_plotAllDirectionsAllCells.m

%% get neuropil signal

load('mask&TCDir.mat');
clear dataTC

npilMask = imCellNeuropil(mask_cell,3,5);
npilTC = stackGetTimeCourses(data_reg,npilMask);

save('npilMask&TCDir.mat','npilMask','npilTC');
%% dF/F
stimOFF_ind = 1:nOFF+nON:size(npilTC,1);

dF_npil = zeros(size(npilTC));
dFoverF_npil = zeros(size(npilTC));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_npil(indAll,:) = bsxfun(@minus,npilTC(indAll,:),mean(npilTC(indF,:),1));
    dFoverF_npil(indAll,:) = bsxfun(@rdivide,dF_npil(indAll,:),mean(npilTC(indF,:),1));
end
%% dF/F (by cell) for each stimulus type

% find on indices for the first frame of each stimulus start period and iti (Off) period

stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_npil,1);

% sort data_TC into 20 frame (10 pre, 10 post) traces around stimON 

dFoverFNpilTrials = zeros(10+nON,size(dFoverF_npil,2),nTrials);
for i = 1:nTrials
    dFoverFNpilTrials(:,:,i) = dFoverF_npil(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

dFoverF_meanDirResp_npil = zeros(size(dFoverFNpilTrials,1),size(dFoverFNpilTrials,2),nStim);
for i = 1:nStim
    trials = find(DirectionDeg == Dirs(i));
    trials = trials(trials<=nTrials);
    dFoverF_meanDirResp_npil(:,:,i) = mean(dFoverFNpilTrials(:,:,trials),3);
end

figure;
for i = 1:nStim
    plot(dFoverF_meanDirResp_npil(:,4,i));
    hold on
end
%% subtract neuropil timecourse from cell timecourse

dFoverF_meanDirResp_sub = dFoverF_meanDirResp - dFoverF_meanDirResp_npil;
dFoverFCellsTrials_sub = dFoverFCellsTrials - dFoverFNpilTrials;

%% plot all cells
cMap = colormap(jet(nStim));
start = 1;
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:15
    cell = start;
    ymax = .1;
    subplot(4,4,iplot);
    for i = 1:nStim
        plot(dFoverF_meanDirResp_sub(:,cell,i),'color',cMap(i,:));
        hold on
        vline(10,'k');
        hold on
        title(['cell' num2str(cell)]);
        hold on
        ymax_i = max(dFoverF_meanDirResp_sub(:,cell,i),[],1);
        if ymax_i > ymax
            ymax = ymax_i;
        end
        axis([0 20 -0.05 ymax]);
        hold on
        legendInfo{i} = [num2str(Dirs(i)) ' degrees'];
        hold on
    end
    if iplot == 1
        legend(legendInfo,'Location','SouthEast')
    end
    start = start+1;
    hold on
end
end

%% plot tuning curve all cells

dFoverFDirResp_sub = zeros(nStim,size(data_TC,2));
errbar = zeros(nStim,size(data_TC,2));
for i = 1:nStim 
    trials = find(DirectionDeg == Dirs(i));
    dFoverFDirResp_sub(i,:) = squeeze(mean(mean(dFoverFCellsTrials_sub(:,:,trials),1),3));
    errbar(i,:) = std(mean(dFoverFCellsTrials_sub(:,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials_sub(:,:,trials),3));
end


start = 1;
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:15
%     cell = baselineStimRespIndex_V(start);
    cell = start;
    ymax = .1;
    subplot(4,4,iplot);
    errorbar(dFoverFDirResp_sub(:,cell),errbar(:,cell),'k');
    hold on
    ymax_i = max(dFoverFDirResp_sub(:,cell),[],1);
    if ymax_i > ymax
        ymax = ymax_i;
    end
    title(['Cell ' num2str(cell)]);
    hold on
    axis([0 9 -0.05 ymax]);
    hold on
%     if ismember(cell,baselineStimRespIndex_V) > 0
%         set(subplot(4,4,iplot),'color',[0.9 0.9 0.9])
%     end
    hold on
    start = start+1;
end
end

