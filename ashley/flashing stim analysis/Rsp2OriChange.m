%% Find frames around orientation change - Fake mouse on
% name and convert some mworks and other variables
cLeverDown = double(cell2mat(input.cLeverDown));
nTrials = (input.trialSinceReset)-1;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
nTrials = input.trialSinceReset-1;
cLeverDown = double(cell2mat(input.cLeverDown));
cTargetOn = input.cTargetOn;
cTargetOn = (double(cell2mat_padded(cTargetOn)))';
cLeverUp = double(cell2mat(input.cLeverUp));
tCyclesOn = double(cell2mat(input.tCyclesOn));
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
Block2ON = Block2ON(:,1:end-1);
TrialOutcome = input.trialOutcomeCell;
Cycles = unique(tCyclesOn);

%pull frames around target stimulus
L = 30;
framesaroundtarget = zeros(size(data_reg,1),size(data_reg,2),nTrials*L);
start = 1;
for itrial = 1:nTrials
    ind = cTargetOn(itrial)-15:cTargetOn(itrial)+14;
    framesaroundtarget(:,:,start:start+L-1) = data_reg(:,:,ind);
    start = start+L;
end
writetiff(framesaroundtarget,'PrePostTarget.tif');

%dF/F and ROI selection
tTrials = size(framesaroundtarget,3)./L;

for itrial = 1:tTrials
    F_byTrial(:,:,itrial) = mean(framesaroundtarget(:,:,(L.*(itrial-1)+1):(L.*(itrial-1)+1)+9),3);
end

dF_data = zeros(size(framesaroundtarget));
for itrial = 1:tTrials
    start1 = L.*(itrial-1)+1;
    start2 = L.*itrial;
    dF_data(:,:,start1:start2) = bsxfun(@minus,framesaroundtarget(:,:,start1:start2),F_byTrial(:,:,itrial));
end

dFoverF_data = zeros(size(framesaroundtarget));
for itrial = 1:tTrials
    start1 = L.*(itrial-1)+1;
    start2 = L.*itrial;
    dFoverF_data(:,:,start1:start2) = bsxfun(@rdivide,dF_data(:,:,start1:start2),F_byTrial(:,:,itrial));
end

last_ind = zeros(1,10);
dF_lastframes_mean_withintrial = zeros(size(dFoverF_data,1),size(dFoverF_data,2),tTrials);
for itrial = 1:tTrials
    last_ind = (L.*itrial)-9:L.*itrial;
    dF_lastframes_mean_withintrial(:,:,itrial) = mean(dFoverF_data(:,:,last_ind),3);
end
clear last_ind

max_dF_lastframes = max(dF_lastframes_mean_withintrial,[],3);
figure; imagesq(max_dF_lastframes); colormap(gray)

bwout = imCellEditInteractive(max_dF_lastframes);
mask_cell = bwlabel(bwout);

data_TC = stackGetTimeCourses(dFoverF_data,mask_cell);
figure; tcOffsetPlot(data_TC)

%find significant responses to vis stim, paired t-test to test visual response
last_ind = zeros(1,10);
data_TC_lastframes = zeros(tTrials,nCells);
for itrial = 1:tTrials
    last_ind = ((L.*itrial)-9:L.*itrial)';
    for icell = 1:nCells
        data_TC_lastframes(itrial,icell) = mean(data_TC(last_ind,icell),1);
    end
end
clear last_ind

first_ind = zeros(1,5);
data_TC_firstframes = zeros(tTrials,nCells);
for itrial = 1:tTrials
    first_ind = (L.*(itrial-1))+10:((L.*(itrial-1))+1)+10;
    for icell = 1:nCells
        data_TC_firstframes(itrial,icell) = mean(data_TC(first_ind,icell),1);
    end
end
clear first_ind

    %paired t-test
first2last_ttest = ttest(data_TC_firstframes,data_TC_lastframes);
visualresp_ind = find(first2last_ttest == 1);

%plot responses and separate by Block2
