SubNum = '516';
date = '141003';
time = '1210';
ImgFolder = '005';
mouse = '516';
fName = '005_000_000';

% load MWorks file
% load MWorks file
CD = ['Z:\data\' mouse '\mworks\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);

fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, ImgFolder);
cd(fileSave);
load('mask&TCDir.mat')
load('oriTuningPreferences.mat')

%%
orig_rate = 30;
final_rate = 3;
down = orig_rate./final_rate;
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
    plot(dFoverF_meanDirResp(:,40,i));
    hold on
end

%% cells' average response to stim
meanCellResp2Stim = zeros(size(dFoverF_meanDirResp,2),nStim);
for i = 1:nStim
    meanCellResp2Stim(:,i) = squeeze(mean(dFoverF_meanDirResp(nOFF+1:end,:,i),1));
end

%% sort cells by preferred direction
meanCellResp2Stim_prefIncluded = meanCellResp2Stim;
meanCellResp2Stim_prefIncluded(:,nStim+1) = dirPref_ind;
meanCellResp2Stim_sortedByPref = sortrows(meanCellResp2Stim_prefIncluded,nStim+1);
meanCellResp2Stim_sortedByPref = meanCellResp2Stim_sortedByPref(:,1:nStim);

%% plot population response
for i = 1:nStim
    totalCellsPref(i) = sum(dirPref_ind == i);
end


figure;
for i = 1:nStim
    subplot(4,3,i)
    plot(meanCellResp2Stim_sortedByPref(:,i));
    hold on
    line = 0;
    xlim([0 size(dFoverF_data,2)]);
    for iline = 1:nStim
        line = line+totalCellsPref(iline);
        vline(line,'k')
        hold on
    end
    dir = Dirs(i);
    stimulus = num2str(dir);
    title([stimulus ' degrees']);
    hold on
    xlabel('cells sorted by direction preference')
    hold on
    ylabel('average dF/F')
    hold on
end


