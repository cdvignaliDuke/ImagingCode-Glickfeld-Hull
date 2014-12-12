data_min = readtiff('minCyclesTrialsPlusAll.tif');
load('tTrialsMinAll.mat');
load('LminPlus.mat');
load('mWorksAll.mat');
load('mask&TCAvg.mat');

data_min = double(data_min);
nframes = tTrialsMinAll*LminPlus;
for i = 1:nframes
    dataSquish(:,:,i) = data_min(:,:,i)*info.S;
end
clear data_min
data_min = dataSquish;
clear dataSquish


for itrial = 1:tTrialsMinAll
    F_byTrial(:,:,itrial) = mean(data_min(:,:,(LminPlus.*(itrial-1)+1):(LminPlus.*(itrial-1)+1)+9),3);
end
dF_data = zeros(size(data_min));
for itrial = 1:tTrialsMinAll
    start1 = LminPlus.*(itrial-1)+1;
    start2 = LminPlus.*itrial;
    dF_data(:,:,start1:start2) = bsxfun(@minus,data_min(:,:,start1:start2),F_byTrial(:,:,itrial));
end

dFoverF_data = zeros(size(data_min));
for itrial = 1:tTrialsMinAll
    start1 = LminPlus.*(itrial-1)+1;
    start2 = LminPlus.*itrial;
    dFoverF_data(:,:,start1:start2) = bsxfun(@rdivide,dF_data(:,:,start1:start2),F_byTrial(:,:,itrial));
end

% dataAvg = readtiff('MinAvgdFoverF.tif');

% bwout = imCellEditInteractive(dataAvg);
% mask_cellAvg = bwlabel(bwout);

data_TCAvg = stackGetTimeCourses(dFoverF_data,mask_cellAvg);
figure; tcOffsetPlot(data_TCAvg)

nCellsAvg = size(data_TCAvg,2);
AvgTrialsCells_mat = zeros(LminPlus,tTrialsMinAll,nCellsAvg);
AvgTrialsCells_trialmean = zeros(LminPlus,nCellsAvg);
for icell = 1:nCellsAvg
    for itrial = 1:tTrialsMinAll
        AvgTrialsCells_mat(:,itrial,icell) = data_TCAvg(1+(LminPlus.*(itrial-1)):LminPlus.*itrial,icell);
    end
    AvgTrialsCells_trialmean(:,icell) = mean(AvgTrialsCells_mat(:,:,icell),2);
end

save('mask&TCAvg.mat','mask_cellAvg','data_TCAvg');

figure;
plot(AvgTrialsCells_trialmean,'k');

AvgTrialsCells_cellsmean = mean(AvgTrialsCells_trialmean,2);

figure; plot(AvgTrialsCells_cellsmean,'k')
axis([0 LminPlus -0.02 0.2]);
ylabel('dF/F');
for i = 1:8
    v = 10+ (i-1)*10.5;
    vline(v,'k');
end

figure;
color = colormap(hsv(nCellsAvg));
for icell = 1:nCellsAvg
    subplot(6,5,icell);
    cellColor = color(icell,:);
    for itrial = 1:tTrialsMinAll
        plot(AvgTrialsCells_mat(:,itrial,icell),'Color',cellColor);
        hold on;
    end
    hold on;
end

A_ind = find(Block2OnSuccess == 1);
V_ind = find(Block2OnSuccess == 0);

AvgTrialsCells_Atrials = squeeze(mean(AvgTrialsCells_mat(:,A_ind,:),2));
AvgTrialsCells_Vtrials = squeeze(mean(AvgTrialsCells_mat(:,V_ind,:),2));

AvgCells_Atrials = squeeze(mean(AvgTrialsCells_Atrials,2));
AvgCells_Vtrials = squeeze(mean(AvgTrialsCells_Vtrials,2));

figure;
plot(AvgCells_Atrials, 'r');
hold on;
plot(AvgCells_Vtrials, 'g');

figure;
color = colormap(hsv(nCellsAvg));
red = color(1,:);
green = color(11,:);
blue = color(19,:);
color = cat(1,red,green,blue);

figure;
start = 1;
for icell = [9 13 15]
    cellColor = color(start,:);
    plot(AvgTrialsCells_trialmean(:,icell),'Color',cellColor);
    hold on;
    axis([0 LminPlus -0.05 0.5]);
    for i = 1:minCyclesOn+1
        v = 10+ (i-1)*10.5;
        vline(v,'k');
    end
    start = start+1;
end

figure;
start = 1;
for icell = [9 13 15]
    cellColor = color(start,:);
    plot(AvgTrialsCells_Vtrials(:,icell),'Color',cellColor);
    hold on;
    axis([0 LminPlus -0.05 0.5]);
    for i = 1:minCyclesOn+1
        v = 10+ (i-1)*10.5;
        vline(v,'k');
    end
    start = start+1;
end

%normalize to peak of first response
AvgTrialsCells_VtrialsNorm = zeros(LminPlus,nCellsAvg);
for icell = 1:nCellsAvg
    ind = max(AvgTrialsCells_Vtrials(10:20,icell));
    AvgTrialsCells_VtrialsNorm(:,icell) = AvgTrialsCells_Vtrials(:,icell)/ind;
end

figure;
start = 1;
for icell = [9 13 15]
    cellColor = color(start,:);
    plot(AvgTrialsCells_VtrialsNorm(:,icell),'Color',cellColor);
    hold on;
%     axis([0 LminPlus -0.05 0.5]);
    for i = 1:minCyclesOn+1
        v = 10+ (i-1)*10.5;
        vline(v,'k');
    end
    start = start+1;
end

figure;
color = colormap(hsv(nCellsAvg));
for icell = 1:nCellsAvg
    subplot(6,5,icell);
    cellColor = color(icell,:);
    for itrial = 1:tTrialsMinAll
        plot(AvgTrialsCells_VtrialsNorm(:,icell),'Color',cellColor);
        hold on;
    end
    hold on;
end