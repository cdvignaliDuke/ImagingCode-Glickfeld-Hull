%% Parameters
SubNum = '510';
date = '140903';
time = '003';
ImgFolder = '003';
mouse = '510';
fName = '003_000_000';
experiment = 'Flashing Stim';


%%

%load retinotopy info
retfile = ['Z:\2P imaging\Analysis\' mouse '\' experiment '\' date '\Retinotopy_V1\analysis']; 
load (retfile);
clear data_reg

%load dataset and mworks file with 'Load_SBXdataset_fast.m'
edit Load_SBXdataset_fast.m

%remove negative data (by addition)
data_sub = data-min(min(min(data,[],1),[],2),[],3);
data = data_sub;
clear data_sub

[out data_reg] = stackRegister(data, data_avg);
clear data

writetiff(data_reg(:,:,1:2000),'FStiff4500.tif');
save ('analysis','-v7.3')
%% call some mworks variables

nTrials = input.trialSinceReset;
cLeverDown = double(cell2mat(input.cLeverDown));
cTargetOn = double(cell2mat(input.cTargetOn));
cLeverUp = double(cell2mat(input.cLeverUp));
tCyclesOn = double(cell2mat(input.tCyclesOn));
ONms = input.stimOnTimeMs;
OFFms = input.stimOffTimeMs;
RateFRperMS = 30./1000;
Block2ON = double(cell2mat(input.tBlock2TrialNumber));
TrialOutcome = input.trialOutcomeCell;
Cycles = unique(tCyclesOn);

%% Find off and on indices for each trial

start = 1;
start1 = 1;
for itrial = 1:nTrials
    start2 = cLeverDown(itrial)-start +1;
    nOFF_ind (1,start1:start1+start2-1) =  start:cLeverDown(itrial);
    start1 = start1+start2;
    start = cLeverUp(itrial);
    
end
siz = size(nOFF_ind,2);
nOFF_ind = nOFF_ind(:,1:siz-2);

%F=average fluorescence for all iti
F1 = mean(data_reg(:,:,nOFF_ind),3);

%% Find dFoverF for each trial for F1 iti

dF1_data = bsxfun(@minus,data_reg, F1);
dF1overF1_data = bsxfun(@rdivide, dF1_data, F1);
max_dF1 = max(dF1overF1_data,[],3);
figure; imagesq(max_dF1); colormap(gray)   

%%
%using Retinotopy_V1 ROI mask, find ROIs in flashingstim dataset
data_TC_FS = stackGetTimeCourses(dF1overF1_data,mask_cell);
figure; tcOffsetPlot(data_TC_FS)

%save data_reg and ROI timecourses
save ('timecourses.mat');

%%
%plot data_TC_FS for just cells that have the preferred retinotopy
CellPrefTimecourses = figure; tcOffsetPlot(data_TC_FS(:,FSpos_ind))
saveas(CellPrefTimecourses,'CellPrefTimecourses.fig')    
    

%%
%Average like-trials for one cell

uniquecycles = unique(tCyclesOn);
numcycles = size(uniquecycles,2);
maxcycles = max(uniquecycles);
triallength = cLeverUp-cLeverDown;
%convert stim ON/OFF times from ms to frames
frON = RateFRperMS.*([0:ONms+OFFms:350*(numcycles)]);
frOFF = RateFRperMS.*([ONms:ONms+OFFms:350*(numcycles)+100]);

%Cell and type of cycle you are interested in
cycletype = 9;
Cell = FSpos_ind(6);
clear Visdata_TC
clear Auddata_TC


Trials = find(tCyclesOn==Cycles(cycletype));
startV = 1;
startA = 1;

for icycle = 1:numcycles
    z = find(tCyclesOn == uniquecycles(icycle));
    avgL_mat(icycle) = ceil(mean(triallength(z)));
end
L = avgL_mat(cycletype);
for itrial = Trials
    y = (cLeverDown(itrial):(cTargetOn(itrial)+(L-length(cLeverDown(itrial):cTargetOn(itrial)))))';
    if Block2ON(itrial) == 1
        Auddata_TC(:,startA) = data_TC_FS(y,Cell);
        startA = startA+1;    
    elseif Block2ON(itrial) == 0
        Visdata_TC(:,startV) = data_TC_FS(y,Cell);
        startV = startV+1;
    end
end

Visdata_TC_mean = mean(Visdata_TC,2);
Auddata_TC_mean = mean(Auddata_TC,2);

%plot auditory vs. visual trials timecourses

figure;
hold on;
subplot(2,1,1);
plot(Visdata_TC_mean);
x = cLeverDown(itrial):cLeverUp (itrial);
Target = find(x == cTargetOn(itrial));
vline(Target,'g')
    for icycle = 1:Cycles(cycletype)
        vline(frON(icycle),'m');
    end
    for icycle = 1:Cycles(cycletype)
        vline(frOFF(icycle),'k');
    end 
title('Visual Trials');
hold on;
subplot(2,1,2);
plot(Auddata_TC_mean);
vline(Target,'g')
for icycle = 1:Cycles(cycletype)
    vline(frON(icycle),'m');
end
for icycle = 1:Cycles(cycletype)
    vline(frOFF(icycle),'k');
end 
title('Auditory Trials');

%%

%find mean pixel strength for each frame in dataset
dF1overF1all_mean = squeeze(mean((mean(dF1overF1_data(:,:,:),1)),2));
%Plot timecourse of full-field for several trials
cycletype = 9;
Cell = FSpos_ind(6);
Trials = find(tCyclesOn==Cycles(cycletype));

FieldofViewTimecourses10 = figure;
Title = ['Field of View Timecourses - ' num2str(Cycles(cycletype)) ' Cycles Before Target Stim'];
title(Title);
hold on
siz1 = size(Trials,2);
start = 1;
for itrial = Trials
    subplot(siz1,1,start);
    y = (cLeverDown(itrial):(cTargetOn(itrial)+(L-length(cLeverDown(itrial):cTargetOn(itrial)))))';
    Target = find(y == cTargetOn(itrial));   
    if Block2ON(itrial) == 1
        plot(dF1overF1all_mean(y));
        title('Auditory')
        vline(Target,'g')
    for icycle = 1:Cycles(cycletype)
        vline(frON(icycle),'m');
    end
    for icycle = 1:Cycles(cycletype)
        vline(frOFF(icycle),'k');
    end
    elseif Block2ON(itrial) == 0
        plot(dF1overF1all_mean(y));
        title('Visual')
        vline(Target,'g')
    for icycle = 1:Cycles(cycletype)
        vline(frON(icycle),'m');
    end
    for icycle = 1:Cycles(cycletype)
        vline(frOFF(icycle),'k');
    end
    end
    start = start+1;
end


saveas(FieldofViewTimecourses10, 'FieldofViewTimecourses10.fig')
%%
%Plot timecourse of one cell for several trials
Trials = [7 8 9 10 11];
siz1 = size(Trials,2);
Cell = 1;
figure;
hold on
start = 1;
for itrial = Trials
    subplot(siz1,1,start);
    x = cLeverDown(itrial):cLeverUp (itrial);
    plot(data_TC_FS(x,FSpos_ind(1,Cell)));
    Target = 22%find(x == cTargetOn(itrial));
    vline(Target,'g')
    for icycle = 1:tCyclesOn(Trial)
        vline(frON(icycle),'m');
    end
    for icycle = 1:tCyclesOn(Trial)
        vline(frOFF(icycle),'k');
    end
    start = start+1;
end

figure
itrial = 1
subplot(2,1,1);x = cLeverDown(itrial):cLeverUp (itrial);
tcOffsetPlot(data_TC_FS(x,FSpos_ind(1,Cell)));
    
    Target = find(x == cTargetOn(itrial));
    vline(Target,'g')
    for icycle = 1:tCyclesOn(Trial)
        vline(frON(icycle),'m');
    end
    for icycle = 1:tCyclesOn(Trial)
        vline(frOFF(icycle),'k');
    end 
    
itrial = 47
subplot(2,1,2);x = cLeverDown(itrial):cLeverUp (itrial);
tcOffsetPlot(data_TC_FS(x,FSpos_ind(1,Cell)));
    
    Target = find(x == cTargetOn(itrial));
    vline(Target,'g')
    for icycle = 1:tCyclesOn(Trial)
        vline(frON(icycle),'m');
    end
    for icycle = 1:tCyclesOn(Trial)
        vline(frOFF(icycle),'k');
    end 
%%
%Plot timecourse of one trial for several cells
siz = size(FSpos_ind,2);
Trial = 10;
figure;
x = cLeverDown(Trial):cLeverUp(Trial);
for icell = 1:siz
    subplot(siz,1,icell);
    plot(data_TC_FS(x,FSpos_ind(1,icell)));
    hold on; 
    for icycle = 1:tCyclesOn(Trial)
        vline(frON(icycle),'m');
    end
    for icycle = 1:tCyclesOn(Trial)
        vline(frOFF(icycle),'k');
    end
end

%%
%save variables for later analysis
save('analysis.mat')