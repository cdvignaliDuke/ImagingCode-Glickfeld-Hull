%save
Save = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis\Rsp2OriChange'];
cd(Save)

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
SuccessTrials_log = strcmp(TrialOutcome,'success');
SuccessTrials_ind = find(SuccessTrials_log == 1);
IgnoreTrials_log = strcmp(TrialOutcome,'ignore');
IgnoreTrials_ind = find(IgnoreTrials_log ==1);
SuccessANDIgnoreTrials_ind = sort([SuccessTrials_ind IgnoreTrials_ind]);


%pull frames around target stimulus - success only
L = 30;
framesaroundtarget = zeros(size(data_reg,1),size(data_reg,2),nTrials*L);
start = 1;
for itrial = SuccessTrials_ind
    ind = cTargetOn(itrial)-15:cTargetOn(itrial)+14;
    framesaroundtarget(:,:,start:start+L-1) = data_reg(:,:,ind);
    start = start+L;
end
writetiff(framesaroundtarget,'PrePostTargetSuccess.tif');

%dF/F and ROI selection
tTrials = size(framesaroundtarget,3)./L;

%find F per trial (need to calculate from pre-trial start data)
% set pre-trial data equal to x
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

data_avg = mean(F_byTrial,3);
figure; imagesq(data_avg); colormap(gray)

max_dF_lastframes = max(dF_lastframes_mean_withintrial,[],3);
figure; imagesq(max_dF_lastframes); colormap(gray)

bwout = imCellEditInteractive(max_dF);
mask_cell = bwlabel(bwout);

data_TC = stackGetTimeCourses(dFoverF_data,mask_cell);
figure; tcOffsetPlot(data_TC)

%with tuning info
data_TC_DirTuningMask_OriChange = stackGetTimeCourses(dFoverF_data,mask_cell_DirTuning);
figure; tcOffsetPlot(data_TC_DirTuningMask_OriChange)

%find significant responses to vis stim, paired t-test to test visual response


nCells = size(data_TC,2);
last_ind = zeros(1,10);
data_TC_lastframes = zeros(size(V_ind,2),nCells);
for itrial = V_ind
    last_ind = ((L.*itrial)-9:L.*itrial)';
    for icell = 1:nCells
        data_TC_lastframes(itrial,icell) = mean(data_TC(last_ind,icell),1);
    end
end
clear last_ind

first_ind = zeros(1,5);
data_TC_firstframes = zeros(size(V_ind,2),nCells);
for itrial = V_ind
    first_ind = (L.*(itrial-1))+11:((L.*(itrial-1))+1)+14;
    for icell = 1:nCells
        data_TC_firstframes(itrial,icell) = mean(data_TC(first_ind,icell),1);
    end
end
clear first_ind

    %paired t-test
first2last_ttest = ttest(data_TC_firstframes,data_TC_lastframes);
visualresp_ind = find(first2last_ttest == 1);

%plot responses and separate by Block2
nCells = size(data_TC,2);
Cell_framesaroundtarget_mat = zeros(L,tTrials,nCells);
Cell_framesaroundtarget_trialmean = zeros(L,nCells);
for icell = 1:nCells
    for itrial = 1:tTrials
        Cell_framesaroundtarget_mat(:,itrial,icell) = data_TC(1+(L.*(itrial-1)):L.*itrial,icell);
    end
    Cell_framesaroundtarget_trialmean(:,icell) = mean(Cell_framesaroundtarget_mat(:,:,icell),2);
end
clear Cell_preLeverDownplus3visstim_mat

figure; plot(Cell_framesaroundtarget_trialmean)
figure; plot(Cell_framesaroundtarget_trialmean(:,visualresp_ind))

%Auditory vs. visual
A_ind = find(Block2ON == 1);
V_ind = find(Block2ON == 0);

Cell_framesaroundtarget_mat_A = Cell_framesaroundtarget_mat(:,A_ind,:);
for icell = 1:nCells
    Cell_framesaroundtarget_mean_A(:,icell) = mean(Cell_framesaroundtarget_mat_A(:,:,icell),2);
end
Cell_framesaroundtarget_mat_V = Cell_framesaroundtarget_mat(:,V_ind,:);
for icell = 1:nCells
    Cell_framesaroundtarget_mean_V(:,icell) = mean(Cell_framesaroundtarget_mat_V(:,:,icell),2);
end

figure; plot(Cell_framesaroundtarget_mean_A(:,nCells),'r')
hold on; plot(Cell_framesaroundtarget_mean_V(:,nCells), 'g')

figure;
start = 1;
for icell = visualresp_ind
    subplot(2,3,start);
    plot(Cell_framesaroundtarget_mean_A(:,icell), 'r');
    hold on;
    plot(Cell_framesaroundtarget_mean_V(:,icell), 'g');
    axis([0 30 -0.1 0.5]);
    vline(15,'k-');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

AllCells_framesaroundtarget_mean_A = mean(Cell_framesaroundtarget_mean_A,2);
AllCells_framesaroundtarget_mean_V = mean(Cell_framesaroundtarget_mean_V,2);

figure; plot(AllCells_framesaroundtarget_mean_A,'r')
hold on; plot(AllCells_framesaroundtarget_mean_V, 'g')
axis([0 30 -0.02 0.3]);
ylabel('dF/F');
vline(15,'k-');
title(['Average All Cells - Auditory vs. Visual']);



%visual trial information
tDirectionDeg = floor(double(cell2mat(input.tGratingDirectionDeg)));

tDirectionDeg_V90_ind = find(tDirectionDeg == 90);
mat_V90 = Cell_framesaroundtarget_mat(:,tDirectionDeg_V90_ind,:);
for icell = 1:nCells
    mean_V90(:,icell) = mean(mat_V90(:,:,icell),2);
end
AllCells_V90 = mean(mean_V90,2);

tDirectionDeg_V63_ind = find(tDirectionDeg == 63);
mat_V63 = Cell_framesaroundtarget_mat(:,tDirectionDeg_V63_ind,:);
for icell = 1:nCells
    mean_V63(:,icell) = mean(mat_V63(:,:,icell),2);
end
AllCells_V63 = mean(mean_V63,2);

tDirectionDeg_V45_ind = find(tDirectionDeg == 45);
mat_V45 = Cell_framesaroundtarget_mat(:,tDirectionDeg_V45_ind,:);
for icell = 1:nCells
    mean_V45(:,icell) = mean(mat_V45(:,:,icell),2);
end
AllCells_V45 = mean(mean_V45,2);

tDirectionDeg_V31_ind = find(tDirectionDeg == 31);
mat_V31 = Cell_framesaroundtarget_mat(:,tDirectionDeg_V31_ind ,:);
for icell = 1:nCells
    mean_V31(:,icell) = mean(mat_V31(:,:,icell),2);
end
AllCells_V31 = mean(mean_V31,2);

tDirectionDeg_V22_ind = find(tDirectionDeg == 22);
mat_V22 = Cell_framesaroundtarget_mat(:,tDirectionDeg_V22_ind,:);
for icell = 1:nCells
    mean_V22(:,icell) = mean(mat_V22(:,:,icell),2);
end
AllCells_V22 = mean(mean_V22,2);

tDirectionDeg_V15_ind = find(tDirectionDeg == 15);
mat_V15 = Cell_framesaroundtarget_mat(:,tDirectionDeg_V15_ind,:);
for icell = 1:nCells
    mean_V15(:,icell) = mean(mat_V15(:,:,icell),2);
end
AllCells_V15 = mean(mean_V15,2);


figure;
plot(mean_V90(:,visualresp_ind),'r');hold on
plot(mean_V63(:,visualresp_ind),'y');hold on
plot(mean_V45(:,visualresp_ind),'g');hold on
plot(mean_V31(:,visualresp_ind),'b');hold on
plot(mean_V22(:,visualresp_ind),'k');hold onvisut
plot(mean_V15(:,visualresp_ind),'m');

% Trial_framesaroundtarget_mat_V63 = Cell_framesaroundtarget_mat(:,tDirectionDeg_above63_ind,:);
% for itrial = tDirectionDeg_above63_ind
%     Cell_framesaroundtarget_mean_V(:,itrial) = mean(Cell_framesaroundtarget_mat(:,itrial,:),3);
% end

Trial_framesaroundtarget_mat_V63 = Cell_framesaroundtarget_mat(:,tDirectionDeg_above63_ind,:);
for itrial = tDirectionDeg_above63_ind
    Cell_framesaroundtarget_mean_Vlessthan(:,itrial) = mean(Cell_framesaroundtarget_mat(:,itrial,:),3);
end

tDirectionDeg_above63_ind = find(tDirectionDeg < 63);

%all trials plotted, all cells averaged.
figure;
start = 1;
for itrial = tDirectionDeg_above63_ind
    subplot(5,5,start);
    plot(Cell_framesaroundtarget_mean_V(:,itrial), 'g');
    axis([0 30 -0.2 0.2]);
    vline(15,'k-');
    title(['Trial ' num2str(itrial)]);
    start = start+1;
end

%one trial, all cells
T = 59;
for icell = 1:nCells
    Trialofinterest(:,icell) = Cell_framesaroundtarget_mat(:,T,icell);
end

figure;
start = 1;
for icell = 1:nCells
    subplot(5,5,start);
    plot(Trialofinterest(:,icell), 'g');
    axis([0 30 -0.2 0.2]);
    vline(15,'k-');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end


%% ******Direction tuning info********
%find significant responses to vis stim, paired t-test to test visual response


nCells = size(data_TC_DirTuningMask_OriChange,2);
last_ind = zeros(1,10);
data_TC_lastframes = zeros(size(V_ind,2),nCells);
for itrial = V_ind
    last_ind = ((L.*itrial)-9:L.*itrial)';
    for icell = 1:nCells
        data_TC_lastframes(itrial,icell) = mean(data_TC_DirTuningMask_OriChange(last_ind,icell),1);
    end
end
clear last_ind

first_ind = zeros(1,5);
data_TC_firstframes = zeros(size(V_ind,2),nCells);
for itrial = V_ind
    first_ind = (L.*(itrial-1))+11:((L.*(itrial-1))+1)+14;
    for icell = 1:nCells
        data_TC_firstframes(itrial,icell) = mean(data_TC_DirTuningMask_OriChange(first_ind,icell),1);
    end
end
clear first_ind

    %paired t-test
first2last_ttest = ttest(data_TC_firstframes,data_TC_lastframes);
visualresp_ind = find(first2last_ttest == 1);

%plot responses and separate by Block2
Cell_framesaroundtarget_mat = zeros(L,tTrials,nCells);
Cell_framesaroundtarget_trialmean = zeros(L,nCells);
for icell = 1:nCells
    for itrial = 1:tTrials
        Cell_framesaroundtarget_mat(:,itrial,icell) = data_TC_DirTuningMask_OriChange(1+(L.*(itrial-1)):L.*itrial,icell);
    end
    Cell_framesaroundtarget_trialmean(:,icell) = mean(Cell_framesaroundtarget_mat(:,:,icell),2);
end
clear Cell_preLeverDownplus3visstim_mat

figure; plot(Cell_framesaroundtarget_trialmean)
figure; plot(Cell_framesaroundtarget_trialmean(:,visualresp_ind))

%Auditory vs. visual
A_ind = find(Block2ON == 1);
V_ind = find(Block2ON == 0);

Cell_framesaroundtarget_mat_A = Cell_framesaroundtarget_mat(:,A_ind,:);
for icell = 1:nCells
    Cell_framesaroundtarget_mean_A(:,icell) = mean(Cell_framesaroundtarget_mat_A(:,:,icell),2);
end
Cell_framesaroundtarget_mat_V = Cell_framesaroundtarget_mat(:,V_ind,:);
for icell = 1:nCells
    Cell_framesaroundtarget_mean_V(:,icell) = mean(Cell_framesaroundtarget_mat_V(:,:,icell),2);
end

figure; plot(Cell_framesaroundtarget_mean_A(:,nCells),'r')
hold on; plot(Cell_framesaroundtarget_mean_V(:,nCells), 'g')

figure;
start = 1;
for icell = [15 19 32]
    subplot(2,3,start);
    plot(Cell_framesaroundtarget_mean_A(:,icell), 'r');
    hold on;
    plot(Cell_framesaroundtarget_mean_V(:,icell), 'g');
    axis([0 30 -0.1 0.5]);
    vline(15,'k-');
    title(['Cell ' num2str(icell)]);
    start = start+1;
end

AllCells_framesaroundtarget_mean_A = mean(Cell_framesaroundtarget_mean_A,2);
AllCells_framesaroundtarget_mean_V = mean(Cell_framesaroundtarget_mean_V,2);

figure; plot(AllCells_framesaroundtarget_mean_A,'r')
hold on; plot(AllCells_framesaroundtarget_mean_V, 'g')
axis([0 30 -0.02 0.3]);
ylabel('dF/F');
vline(15,'k-');
title(['Average All Cells - Auditory vs. Visual']);
