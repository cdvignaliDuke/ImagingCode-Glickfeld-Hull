date = '141003';
mouse = '516';
ImgFolder = '002+003+004';
CD = ['Z:\2P imaging\Analysis\' mouse '\' date '\' ImgFolder '\FlashingStimAnalysis'];
cd(CD);

data_min = readtiff('minCyclesTrialsAll.tif');
data_max = readtiff('maxCyclesTrialsAll.tif');
load('tTrialsMinAll.mat')
load('tTrialsMaxAll.mat')
load('mWorksAll.mat')
load('Lmin.mat')
load('Lmax.mat')

load('mask&TCMin.mat')
load('mask&TCMax.mat')

%% dF/F 

%min

for itrial = 1:tTrialsMinAll
    F_byTrial(:,:,itrial) = mean(data_min(:,:,(Lmin.*(itrial-1)+1):(Lmin.*(itrial-1)+1)+9),3);
end
dF_data = zeros(size(data_min));
for itrial = 1:tTrialsMinAll
    start1 = Lmin.*(itrial-1)+1;
    start2 = Lmin.*itrial;
    dF_data(:,:,start1:start2) = bsxfun(@minus,data_min(:,:,start1:start2),F_byTrial(:,:,itrial));
end

dFoverF_data = zeros(size(data_min));
for itrial = 1:tTrialsMinAll
    start1 = Lmin.*(itrial-1)+1;
    start2 = Lmin.*itrial;
    dFoverF_data(:,:,start1:start2) = bsxfun(@rdivide,dF_data(:,:,start1:start2),F_byTrial(:,:,itrial));
end
    
%get average image for response to vis stim

dF_meanVisRspMin = zeros(size(dFoverF_data,1),size(dFoverF_data,2),tTrialsMinAll);
for itrial = 1:tTrialsMinAll
    last_ind = 11+(Lmin*(itrial-1)):Lmin.*itrial;
    dF_meanVisRspMin(:,:,itrial) = mean(dFoverF_data(:,:,last_ind),3);
end
clear last_ind

max_dF_VisRspMin = max(dF_meanVisRspMin,[],3);
figure; imagesq(max_dF_VisRspMin); colormap(gray)

bwout = imCellEditInteractive(max_dF_VisRspMin);
mask_cellMin = bwlabel(bwout);

data_TCMin = stackGetTimeCourses(dFoverF_data,mask_cellMin);
figure; tcOffsetPlot(data_TCMin)

nCellsMin = size(data_TCMin,2);
MinTrialsCells_mat = zeros(Lmin,tTrialsMinAll,nCellsMin);
MinTrialsCells_trialmean = zeros(Lmin,nCellsMin);
for icell = 1:nCellsMin
    for itrial = 1:tTrialsMinAll
        MinTrialsCells_mat(:,itrial,icell) = data_TCMin(1+(Lmin.*(itrial-1)):Lmin.*itrial,icell);
    end
    MinTrialsCells_trialmean(:,icell) = mean(MinTrialsCells_mat(:,:,icell),2);
end

save('mask&TCMin.mat','mask_cellMin','data_TCMin');
 
%max
for itrial = 1:tTrialsMaxAll
    F_byTrial(:,:,itrial) = mean(data_max(:,:,(Lmax.*(itrial-1)+1):(Lmax.*(itrial-1)+1)+9),3);
end
dF_data = zeros(size(data_max));
for itrial = 1:tTrialsMaxAll
    start1 = Lmax.*(itrial-1)+1;
    start2 = Lmax.*itrial;
    dF_data(:,:,start1:start2) = bsxfun(@minus,data_max(:,:,start1:start2),F_byTrial(:,:,itrial));
end

dFoverF_data = zeros(size(data_max));
for itrial = 1:tTrialsMaxAll
    start1 = Lmax.*(itrial-1)+1;
    start2 = Lmax.*itrial;
    dFoverF_data(:,:,start1:start2) = bsxfun(@rdivide,dF_data(:,:,start1:start2),F_byTrial(:,:,itrial));
end
    
%get average image for response to vis stim

dF_meanVisRspMax = zeros(size(dFoverF_data,1),size(dFoverF_data,2),tTrialsMaxAll);
for itrial = 1:tTrialsMaxAll
    last_ind = 11+(Lmax*(itrial-1)):Lmax.*itrial;
    dF_meanVisRspMax(:,:,itrial) = mean(dFoverF_data(:,:,last_ind),3);
end
clear last_ind

max_dF_VisRspMax = max(dF_meanVisRspMax,[],3);
figure; imagesq(max_dF_VisRspMax); colormap(gray)

bwout = imCellEditInteractive(max_dF_VisRspMax);
mask_cellMax = bwlabel(bwout);

data_TCMax = stackGetTimeCourses(dFoverF_data,mask_cellMax);
figure; tcOffsetPlot(data_TCMax)

nCellsMax = size(data_TCMax,2);
MaxTrialsCells_mat = zeros(Lmax,tTrialsMaxAll,nCellsMax);
MaxTrialsCells_trialmean = zeros(Lmax,nCellsMax);
for icell = 1:nCellsMax
    for itrial = 1:tTrialsMaxAll
        MaxTrialsCells_mat(:,itrial,icell) = data_TCMax(1+(Lmax.*(itrial-1)):Lmax.*itrial,icell);
    end
    MaxTrialsCells_trialmean(:,icell) = mean(MaxTrialsCells_mat(:,:,icell),2);
end

save('mask&TCMax.mat','mask_cellMax','data_TCMax');

%% ** visualizing**
figure;
plot(MinTrialsCells_trialmean,'k');
hold on;
plot(MaxTrialsCells_trialmean,'b');

MinTrialsCells_cellsmean = mean(MinTrialsCells_trialmean,2);
MaxTrialsCells_cellsmean = mean(MaxTrialsCells_trialmean,2);

figure; plot(MinTrialsCells_cellsmean,'k')
hold on; plot(MaxTrialsCells_cellsmean, 'b')
axis([0 Lmax -0.02 0.2]);
ylabel('dF/F');
for i = 1:8
    v = 10+ (i-1)*10.5;
    vline(v,'k');
end
title(['Average All Cells - Long vs. Short']);

figure;
color = colormap(hsv(nCellsMin));
for icell = 1:nCellsMin
    subplot(6,5,icell);
    cellColor = color(icell,:);
    for itrial = 1:tTrialsMinAll
        plot(MinTrialsCells_mat(:,itrial,icell),'Color',cellColor);
        hold on;
    end
    hold on;
end

figure;
color = colormap(hsv(nCellsMin));
for icell = 1:nCellsMin
    cellColor = color(icell,:);
    plot(MinTrialsCells_trialmean(:,icell),'Color',cellColor);
    hold on;
end

% min trials A vs V
minCyclesOn = inputAll(1).minCyclesOn;
start = 1;
for i = 1:3
    siz = size(inputAll(i).trialOutcomeCell,2);
    TrialOutcome(1,start:start+siz-1) = inputAll(i).trialOutcomeCell;
    start = start+siz;
end
success_log = strcmp(TrialOutcome,'success');
success_ind = find(success_log == 1);

start = 1;
for i = 1:3
    siz = size(inputAll(i).tCyclesOn,2);
    tCyclesOn(1,start:start+siz-1) = inputAll(i).tCyclesOn;
    start = start+siz;
end
tCyclesOn = cell2mat(tCyclesOn);
minCycles_ind = find(tCyclesOn == minCyclesOn);
minCycles_trials = intersect(success_ind,minCycles_ind);

start = 1;
for i = 1:3
    siz = size(inputAll(i).tBlock2TrialNumber,2);
    Block2On(1,start:start+siz-1) = inputAll(i).tBlock2TrialNumber;
    start = start+siz;
end
Block2On = cell2mat(Block2On);
Block2OnSuccess = Block2On(:,minCycles_trials);

A_ind = find(Block2OnSuccess == 1);
V_ind = find(Block2OnSuccess == 0);

MinTrialsCells_Atrials = squeeze(mean(MinTrialsCells_mat(:,A_ind,:),2));
MinTrialsCells_Vtrials = squeeze(mean(MinTrialsCells_mat(:,V_ind,:),2));

figure;
for icell = 1:nCellsMin
    subplot(3,4,icell);
    cellColor = color(icell,:);
    plot(MinTrialsCells_Atrials(:,icell),'r');
    hold on;
    plot(MinTrialsCells_Vtrials(:,icell),'g');
    hold on;
    axis([0 Lmin -0.05 0.35]);
    for i = 1:minCyclesOn+1
        v = 10+ (i-1)*10.5;
        vline(v,'k');
    end
end
