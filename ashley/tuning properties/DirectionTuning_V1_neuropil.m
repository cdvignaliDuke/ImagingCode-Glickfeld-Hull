%% get dF/F of cells from direction tuning
edit DirectionTuning_V1_plotAllDirectionsAllCells.m

%% get neuropil signal

load('mask&TCDir.mat');
clear dataTC data_TC

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

%% subtract neuropil from timecourse
dataTC_sub = data_TC - npilTC;
dataTCsubMin = min(min(dataTC_sub,[],2),[],1);
dataTC_sub = dataTC_sub - dataTCsubMin;

%% dF/F (by cell) of npil for each stimulus type

% find on indices for the first frame of each stimulus start period and iti (Off) period

stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_npil,1);

% sort npil into 20 frame (10 pre, 10 post) traces around stimON 

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

%% dF/F (by cell) of dataTC_sub for each stimulus type

stimOFF_ind = 1:nOFF+nON:size(dataTC_sub,1);

dF_dataTCsub = zeros(size(dataTC_sub));
dFoverF_dataTCsub = zeros(size(dataTC_sub));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_dataTCsub(indAll,:) = bsxfun(@minus,dataTC_sub(indAll,:),mean(dataTC_sub(indF,:),1));
    dFoverF_dataTCsub(indAll,:) = bsxfun(@rdivide,dF_dataTCsub(indAll,:),mean(dataTC_sub(indF,:),1));
end

% find on indices for the first frame of each stimulus start period and iti (Off) period

stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_dataTCsub,1);

% sort npil into 20 frame (10 pre, 10 post) traces around stimON 

dFoverFdataTCsubTrials = zeros(10+nON,size(dFoverF_dataTCsub,2),nTrials);
for i = 1:nTrials
    dFoverFdataTCsubTrials(:,:,i) = dFoverF_dataTCsub(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

dFoverF_meanDirResp_dataTCsub = zeros(size(dFoverFdataTCsubTrials,1),size(dFoverFdataTCsubTrials,2),nStim);
for i = 1:nStim
    trials = find(DirectionDeg == Dirs(i));
    trials = trials(trials<=nTrials);
    dFoverF_meanDirResp_dataTCsub(:,:,i) = mean(dFoverFdataTCsubTrials(:,:,trials),3);
end

figure;
for i = 1:nStim
    plot(dFoverF_meanDirResp_dataTCsub(:,4,i));
    hold on
end

%% plot all neuropil
cMap = colormap(jet(nStim));
start = 1;
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:15
    cell = start;
    ymax = .1;
    subplot(4,4,iplot);
    for i = 1:nStim
        plot(dFoverF_meanDirResp_npil(:,cell,i),'color',cMap(i,:));
        hold on
        vline(10,'k');
        hold on
        title(['neuropil' num2str(cell)]);
        hold on
        ymax_i = max(dFoverF_meanDirResp_npil(:,cell,i),[],1);
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
        plot(dFoverF_meanDirResp_dataTCsub(:,cell,i),'color',cMap(i,:));
        hold on
        vline(10,'k');
        hold on
        title(['cell - npil' num2str(cell)]);
        hold on
        ymax_i = max(dFoverF_meanDirResp_dataTCsub(:,cell,i),[],1);
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

dFoverFDirResp_dataTCsub = zeros(nStim,size(data_TC,2));
errbar = zeros(nStim,size(data_TC,2));
for i = 1:nStim 
    trials = find(DirectionDeg == Dirs(i));
    dFoverFDirResp_dataTCsub(i,:) = squeeze(mean(mean(dFoverFdataTCsubTrials(:,:,trials),1),3));
    errbar(i,:) = std(mean(dFoverFdataTCsubTrials(:,:,trials),1),[],3)/sqrt(size(dFoverFdataTCsubTrials(:,:,trials),3));
end


start = 1;
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:15
%     cell = baselineStimRespIndex_V(start);
    cell = start;
    ymax = .1;
    subplot(4,4,iplot);
    errorbar(dFoverFDirResp_dataTCsub(:,cell),errbar(:,cell),'k');
    hold on
    ymax_i = max(dFoverFDirResp_dataTCsub(:,cell),[],1);
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

