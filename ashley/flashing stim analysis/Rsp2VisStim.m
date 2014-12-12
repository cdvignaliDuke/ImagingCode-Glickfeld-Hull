%load dataset
% SubNum = '510';
% date = '140909';
% time = '1147';
% ImgFolder = '008';
% mouse = '510';
% fName = '008_000_000';
% experiment = 'Flashing Stim';
% 
% load MWorks file
CD = ['Z:\2P imaging\MWorks\' mouse '\' date];
cd(CD);
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (mworks);


%set current directory for saving
%create folders
CD = ['Z:\2P imaging\Analysis\' mouse '\' date ];
cd(CD);
mkdir(ImgFolder);
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder];
cd(CD);
mkdir('FlashingStimAnalysis','Rsp2VisStim');
%save path
Save = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2VisStim'];
cd(Save)

%%
%name and convert some mworks variables
cLeverDown = double(cell2mat(input.cLeverDown));
nTrials = (input.trialSinceReset);
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
nTrials = input.trialSinceReset-1;
cLeverDown = double(cell2mat(input.cLeverDown));
cTargetOn = input.cTargetOn;
for itrial = 1:nTrials
    if isempty(cTargetOn{itrial})
        cTargetOn{itrial} = NaN;
    end
end
cTargetOn = (double(cell2mat_padded(cTargetOn)))'; %For now NaNs == 0, may need to change...
cLeverUp = double(cell2mat(input.cLeverUp));
tCyclesOn = double(cell2mat(input.tCyclesOn));
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
Cycles = unique(tCyclesOn);


% % for fake mouse skip last trial
% cLeverDown = double(cell2mat(input.cLeverDown));
% nTrials = (input.trialSinceReset)-1;
% RateFRperMS = 30./1000;
% Block2ON = double(cell2mat(input.tBlock2TrialNumber));
% Block2ON = Block2ON(:,1:end-1);
% cLeverDown = double(cell2mat(input.cLeverDown));
% cTargetOn = input.cTargetOn;
% for itrial = 1:nTrials
%     if isempty(cTargetOn{itrial})
%         cTargetOn{itrial} = NaN;
%     end
% end
% cTargetOn = (double(cell2mat_padded(cTargetOn)))'; %For now NaNs == 0, may need to change...
% cLeverUp = double(cell2mat(input.cLeverUp));
% tCyclesOn = double(cell2mat(input.tCyclesOn));
% ONms = input.stimOnTimeMs;
% OFFms = input.stimOffTimeMs;
% RateFRperMS = 30./1000;
% Block2ON = double(cell2mat(input.tBlock2TrialNumber));
% TrialOutcome = input.trialOutcomeCell;
% Cycles = unique(tCyclesOn);


%% Tif for ImageJ Analysis
%align frames for just those around the start of each trial
L = 30;
framesaroundleverdown = zeros(size(data_reg,1),size(data_reg,2),nTrials*L);
start = 1;
for itrial = 1:nTrials
    ind = cLeverDown(itrial)-15:cLeverDown(itrial)+14;
    framesaroundleverdown(:,:,start:start+L-1) = data_reg(:,:,ind);
    start = start+L;
end
writetiff(framesaroundleverdown,'PrePostTrialStart.tif');

%% Timecourse of above Tif


%find timecourse for iti right before start of trial
meanFrameData = squeeze(mean(mean(data_reg,1),2));

ITItimecourse = zeros(nTrials,15);
for itrial = 1:nTrials
    ind = cLeverDown(itrial)-15:cLeverDown(itrial)-1;
    ITItimecourse(itrial,:) = meanFrameData(ind);
end
ITItimecourse_mean = mean(ITItimecourse,1);

%find timecourse for first 15 frames of trial
TrStarttimecourse = zeros(nTrials,15);
for itrial = 1:nTrials
    ind = cLeverDown(itrial):cLeverDown(itrial)+14;
    TrStarttimecourse(itrial,:) = meanFrameData(ind);
end
TrStarttimecourse_mean = mean(TrStarttimecourse,1);

PreANDPostTrStartTimecourses = figure;
title('Pre (b) and Post (r) Trial Start');
xlabel('raw fluorescence')
ylabel('frame')
plot(ITItimecourse_mean,'b');
hold on
plot(TrStarttimecourse_mean,'r');

%% Tif for just success and ignore trials

% for real behavior
% SuccessTrials_log = strcmp(TrialOutcome,'success');
% SuccessTrials_ind = find(SuccessTrials_log == 1);
% IgnoreTrials_log = strcmp(TrialOutcome,'ignore');
% IgnoreTrials_ind = find(IgnoreTrials_log ==1);
% SuccessANDIgnoreTrials_ind = sort([SuccessTrials_ind IgnoreTrials_ind]);
% 
% % 10 frames pre-LeverDown plus 3 vis stim (equal to about 32 frames) to
% % capture visual reponses before target.
% L = 10+ceil((RateFRperMS*350)*3);
% siz = size(SuccessANDIgnoreTrials_ind,2)-1;
% preLeverDownplus3visstim = zeros(size(data_reg,1),size(data_reg,2),siz*L);
% start = 1;
% for itrial = SuccessANDIgnoreTrials_ind
%     ind = cLeverDown(itrial)-10:cLeverDown(itrial)+(L-11);
%     preLeverDownplus3visstim(:,:,start:start+L-1) = data_reg(:,:,ind);
%     start = start+L;
% end

% writetiff(preLeverDownplus3visstim,'Rsp2VisStim.tif');
% preLDP3VS04 = preLeverDownplus3visstim; 
% clear preLeverDownplus3visstim
% clear data_reg

% for 'Fake Mouse success only'
L = 10+ceil((RateFRperMS*350)*3);
preLeverDownplus3visstim = zeros(size(data_reg,1),size(data_reg,2),nTrials*L);
start = 1;
for itrial = 1:nTrials
    ind = cLeverDown(itrial)-10:cLeverDown(itrial)+(L-11);
    preLeverDownplus3visstim(:,:,start:start+L-1) = data_reg(:,:,ind);
    start = start+L;
end

writetiff(preLeverDownplus3visstim,'3VisStim.tif');


% preLeverDownplus3visstim_2 = preLeverDownplus3visstim_007;
% clear preLeverDownplus3visstim


% preLeverDownplus3visstim_007 = zeros(size(preLeverDownplus3visstim_1,1),size(preLeverDownplus3visstim_1,2),(size(preLeverDownplus3visstim_1,3)+size(preLeverDownplus3visstim_2,3)));
% preLeverDownplus3visstim_007 = cat(3,preLeverDownplus3visstim_1, preLeverDownplus3visstim_2);
% clear preLeverDownplus3visstim_1;
% % clear preLeverDownplus3visstim_2;
% writetiff(preLeverDownplus3visstim_005and006,'preLeverDownplus3visstim_005+006.tif');
%% Cell responses
% siz1 = size(preLeverDownplus3visstim_005,1);
% siz2 = size(preLeverDownplus3visstim_005,2);
tTrials = size(preLeverDownplus3visstim,3)./L;
% F_ind = zeros(1,tTrials.*10);
% start = 1;
% for itrial = 1:tTrials
%     F_ind(start:start+9) = (L.*(itrial-1)+1):(L.*(itrial-1)+1)+9;
%     start = start+10;
% end

% 
% F = mean(preLeverDownplus3visstim_005and006(:,:,F_ind),3);
% dF_data = bsxfun(@minus,preLeverDownplus3visstim_005and006, F);
% dFoverF_data = bsxfun(@rdivide, dF_data, F);

for itrial = 1:tTrials
    F_byTrial(:,:,itrial) = mean(preLeverDownplus3visstim(:,:,(L.*(itrial-1)+1):(L.*(itrial-1)+1)+9),3);
end
dF_data = zeros(size(preLeverDownplus3visstim));
for itrial = 1:tTrials
    start1 = L.*(itrial-1)+1;
    start2 = L.*itrial;
    dF_data(:,:,start1:start2) = bsxfun(@minus,preLeverDownplus3visstim(:,:,start1:start2),F_byTrial(:,:,itrial));
end

dFoverF_data = zeros(size(preLeverDownplus3visstim));
for itrial = 1:tTrials
    start1 = L.*(itrial-1)+1;
    start2 = L.*itrial;
    dFoverF_data(:,:,start1:start2) = bsxfun(@rdivide,dF_data(:,:,start1:start2),F_byTrial(:,:,itrial));
end
    
%get average image for last 20 of all frames(42) per trial
last_ind = zeros(1,20);
dF_lastframes_mean_withintrial = zeros(size(dFoverF_data,1),size(dFoverF_data,2),tTrials);
for itrial = 1:tTrials
    last_ind = (L.*itrial)-19:L.*itrial;
    dF_lastframes_mean_withintrial(:,:,itrial) = mean(dFoverF_data(:,:,last_ind),3);
end
clear last_ind

dF_lastframes = mean(dFoverF_data(:,:,last_ind),3);
figure; imagesq(dF_lastframes); colormap(gray)

max_dF_lastframes = max(dF_lastframes_mean_withintrial,[],3);
figure; imagesq(max_dF_lastframes); colormap(gray)

% max_dF = max(dFoverF_data,[],3);
% figure; imagesq(max_dF); colormap(gray)

bwout = imCellEditInteractive(max_dF_lastframes);
mask_cell = bwlabel(bwout);

%used average F for all trials to get a brighter max_dF for imCellEdit, now
%re-baseline to F on a trial by trial basis

% dF_data = zeros(size(preLeverDownplus3visstim_007));
% for itrial = 1:tTrials
%     start1 = L.*(itrial-1)+1;
%     start2 = L.*itrial;
%     dF_data(:,:,start1:start2) = bsxfun(@minus,preLeverDownplus3visstim_007(:,:,start1:start2),F_byTrial(:,:,itrial));
% end
% 
% dFoverF_data = zeros(size(preLeverDownplus3visstim_007));
% for itrial = 1:tTrials
%     start1 = L.*(itrial-1)+1;
%     start2 = L.*itrial;
%     dFoverF_data(:,:,start1:start2) = bsxfun(@rdivide,dF_data(:,:,start1:start2),F_byTrial(:,:,itrial));
% end

data_TC = stackGetTimeCourses(dFoverF_data,mask_cell);
figure; tcOffsetPlot(data_TC)

data_TC_DirTuningMask = stackGetTimeCourses(dFoverF_data,mask_cell_DirTuning);
figure; tcOffsetPlot(data_TC_DirTuningMask)

%% 

%for one cell, average all trials

% % for real behavior
% % 10 frames pre-LeverDown plus 3 vis stim (equal to about 32 frames) to
% % capture visual reponses before target.
% L = 10+ceil((RateFRperMS*350)*3);
% siz = size(SuccessANDIgnoreTrials_ind,2);
% Cell4_preLeverDownplus3visstim = zeros(1,siz*L);
% Cell = 6;
% start = 1;
% for itrial = SuccessANDIgnoreTrials_ind
%     ind = cLeverDown(itrial)-10:cLeverDown(itrial)+(L-11);
%     Cell4_preLeverDownplus3visstim(start:start+L-1) = data_TC(ind,Cell);
%     start = start+L;
% end

% for 'Fake mouse success only'

nCells = size(data_TC,2);
Cell_preLeverDownplus3visstim_mat = zeros(L,tTrials,nCells);
Cell_preLeverDownplus3visstim_trialmean = zeros(L,nCells);
for icell = 1:nCells
    for itrial = 1:tTrials
        Cell_preLeverDownplus3visstim_mat(:,itrial,icell) = data_TC(1+(42.*(itrial-1)):42.*itrial,icell);
    end
    Cell_preLeverDownplus3visstim_trialmean(:,icell) = mean(Cell_preLeverDownplus3visstim_mat(:,:,icell),2);
end
clear Cell_preLeverDownplus3visstim_mat

figure; plot(Cell_preLeverDownplus3visstim_trialmean)

%find significant responses to vis stim, paired t-test to test visual response
last_ind = zeros(1,20);
data_TC_lastframes = zeros(tTrials,nCells);
for itrial = 1:tTrials
    last_ind = ((L.*itrial)-19:L.*itrial)';
    for icell = 1:nCells
        data_TC_lastframes(itrial,icell) = mean(data_TC(last_ind,icell),1);
    end
end
clear last_ind

first_ind = zeros(1,10);
data_TC_firstframes = zeros(tTrials,nCells);
for itrial = 1:tTrials
    first_ind = (L.*(itrial-1))+1:((L.*(itrial-1))+1)+10;
    for icell = 1:nCells
        data_TC_firstframes(itrial,icell) = mean(data_TC(first_ind,icell),1);
    end
end
clear first_ind

    %paired t-test
first2last_ttest = ttest(data_TC_firstframes,data_TC_lastframes);
visualresp_ind = find(first2last_ttest == 1);
% first2last_ttest = ttest(data_TC_firstframes,data_TC_lastframes);
% visualresp_ind_FbyTrial = find(first2last_ttest == 1);

% figure; scatter(ones(1,tTrials),data_TC_firstframes(:,2),'k'); 
% hold on; scatter((2.*(ones(1,tTrials))),data_TC_lastframes(:,2),'b')

    %plot
figure; plot(Cell_preLeverDownplus3visstim_trialmean(:,visualresp_ind))
% figure; plot(Cell_preLeverDownplus3visstim_trialmean(:,visualresp_ind_FbyTrial))


save('Rsp2VisStimAnalysis.mat', 'mask_cell', 'data_TC', 'visualresp_ind');

%% Direction tuning cell mask
last_ind = zeros(1,20);
data_TC_lastframes = zeros(tTrials,nCells);
for itrial = 1:tTrials
    last_ind = ((L.*itrial)-19:L.*itrial)';
    for icell = 1:nCells
        data_TC_lastframes(itrial,icell) = mean(data_TC_DirTuningMask(last_ind,icell),1);
    end
end
clear last_ind

first_ind = zeros(1,10);
data_TC_firstframes = zeros(tTrials,nCells);
for itrial = 1:tTrials
    first_ind = (L.*(itrial-1))+1:((L.*(itrial-1))+1)+10;
    for icell = 1:nCells
        data_TC_firstframes(itrial,icell) = mean(data_TC_DirTuningMask(first_ind,icell),1);
    end
end
clear first_ind

first2last_ttest = ttest(data_TC_firstframes,data_TC_lastframes);
visualresp_ind = find(first2last_ttest == 1);
