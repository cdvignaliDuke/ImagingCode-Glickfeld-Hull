SubNum = '613';
date = '150511';
time = '1513';
ImgFolder = '005';
mouse = 'AW13';

% load MWorks file
% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

% %load preferences
% fsFolder = '004+005';
% fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, fsFolder);
% cd(fileSave);
% load('cellSelectivityIndices.mat');

% %load data_TC
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
% load('mask&TCDir.mat')
% % load('dataTC.mat')
% edit DirectionTuning_V1.m
load('neuropil.mat')
data_TC = npSubTC;
%%
down = 10;
nON = (input.nScansOn)./down;
nOFF = (input.nScansOff)./down;
nStim = input.gratingDirectionStepN;
nRep = size(data_TC,1)./((nON+nOFF)*nStim);
nTrials = (nStim.*nRep);
DirectionDeg = cell2mat(input.tGratingDirectionDeg);
Dirs = unique(DirectionDeg);

%% dF/F by trial
% F per trial
stimOFF_ind = 1:nOFF+nON:size(data_TC,1);

dF_data = zeros(size(data_TC));
dFoverF_data = zeros(size(data_TC));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_data(indAll,:) = bsxfun(@minus,data_TC(indAll,:),mean(data_TC(indF,:),1));
    dFoverF_data(indAll,:) = bsxfun(@rdivide,dF_data(indAll,:),mean(data_TC(indF,:),1));
end
%% dF/F (by cell) for each stimulus type

% find on indices for the first frame of each stimulus start period and iti (Off) period

stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_data,1);

% sort data_TC into 20 frame (10 pre, 10 post) traces around stimON 

dFoverFCellsTrials = zeros(10+nON,size(dFoverF_data,2),nTrials);
for i = 1:nTrials
    dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

dFoverF_meanDirResp = zeros(size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2),nStim);
for i = 1:nStim
    trials = find(DirectionDeg == Dirs(i));
    dFoverF_meanDirResp(:,:,i) = mean(dFoverFCellsTrials(:,:,trials),3);
end

figure;
for i = 1:nStim
    plot(dFoverF_meanDirResp(:,10,i));
    hold on
end

%% 
cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);

cMap = colormap(jet(nStim));
start = 1;
errbar = zeros(nStim,size(data_TC,2));
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:31
%     cell = baselineStimRespIndex_V(start);
    cell = start;
    ymax = .1;
    subplot(8,4,iplot);
    for i = 1:nStim
        plot(dFoverF_meanDirResp(:,cell,i),'color',cMap(i,:));
        hold on
        vline(10,'k');
        hold on
        title(['Cell ' num2str(cell)]);
        hold on
        ymax_i = max(dFoverF_meanDirResp(:,cell,i),[],1);
        if ymax_i > ymax
            ymax = ymax_i;
        end
        axis([0 20 -0.05 ymax]);
        hold on
        legendInfo{i} = [num2str(Dirs(i)) ' degrees'];
        hold on

    end   
    if ismember(cell,cellsPrefZero) > 0
        set(subplot(8,4,iplot),'color',[0.9 0.9 0.9])
    end
    if ismember(cell,cellsPrefNinety) > 0
        set(subplot(8,4,iplot),'color',[0.8 0.8 0.8])
    end
    if iplot == 15 & ifig == 1
        legend(legendInfo,'Location','SouthEast')
    end
    start = start+1;
    hold on
end
end

dFoverFDirResp = zeros(nStim,size(data_TC,2));
errbar = zeros(nStim,size(data_TC,2));
for i = 1:nStim 
    trials = find(DirectionDeg == Dirs(i));
    dFoverFDirResp(i,:) = squeeze(mean(mean(dFoverFCellsTrials(11:16,:,trials),1),3));
    errbar(i,:) = std(mean(dFoverFCellsTrials(11:16,:,trials),1),[],3)/sqrt(size(dFoverFCellsTrials(11:16,:,trials),3));
end


start = 1;
for ifig = 1:ceil(size(data_TC,2)/15)
figure;
for iplot = 1:31
%     cell = baselineStimRespIndex_V(start);
    cell = start;
    ymax = .1;
    subplot(8,4,iplot);
    errorbar(dFoverFDirResp(:,cell),errbar(:,cell),'k');
    hold on
    ymax_i = max(dFoverFDirResp(:,cell),[],1);
    if ymax_i > ymax
        ymax = ymax_i;
    end
    title(['Cell ' num2str(cell)]);
    hold on
    axis([0 length(Dirs) -0.05 ymax]);
    hold on
%     if ismember(cell,baselineStimRespIndex_V) > 0
%         set(subplot(4,4,iplot),'color',[0.9 0.9 0.9])
%     end
    hold on
    start = start+1;
end
end

%% select cells
cells = [ 9 165];
figure;
for iplot = 1:9
%     cell = baselineStimRespIndex_V(start);
    cell = cells(iplot);
    ymax = .1;
    subplot(3,3,iplot);
    errorbar(dFoverFDirResp(:,cell),errbar(:,cell),'k');
    hold on
    ymax_i = max(dFoverFDirResp(:,cell),[],1);
    if ymax_i > ymax
        ymax = ymax_i;
    end
    title(['Cell ' num2str(cell)]);
    hold on
    axis([0 length(Dirs) -0.05 ymax]);
    hold on
end
%%
cells = [33 121];
col = {'y','c'};
legendInfo = [];
figure;
for icell = 1:length(cells)
    errorbar(Dirs,dFoverFDirResp(:,cells(icell)),errbar(:,cells(icell)),[col{icell} 'o'],'MarkerSize', 10)
    legendInfo{icell} = ['Cell ' num2str(cells(icell))];
    hold on
end
ylim([-0.1 0.4]) 
xlim([-10 Dirs(end)+5])
xlabel('Direction')
ylabel('dF/F')
set(gca,'XTick',Dirs)
legend(legendInfo)
title('Select Cell Tuning')
axis square

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

print(['Z:\Analysis\temp figs\150924retreatposter' '\'  date mouse '_selectcelltuning'], '-dpdf')

%     print([fnout '\'  date '_' figNameBase '_STDrsplastcyc_audVSvis'], '-dpdf'


    
    
    