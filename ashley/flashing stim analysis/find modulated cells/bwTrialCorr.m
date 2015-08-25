run('FlashingStim_dataSortedByCycle_combineDatasets.m')
%%
dirFolder = '006';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, dirFolder);
cd(fileSave);
% load('oriTuningPreferences.mat')
load('TuningPreferences.mat')


%%
% cell types
dataTrialStart = cycDataDFoverF_cmlvNoTarget{4};
v_ind = cycV_ind{4};
% a_ind = cycA_ind{1};
a_ind = cycAV_ind{4};


preStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial =1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        preStimResp_V(itrial,icell) = mean(dataTrialStart(1:30,icell,v_ind(itrial)),1);
    end
end

baselineStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
for itrial = 1:size(v_ind,1);
    for icell = 1:size(dataTrialStart,2)
        baselineStimResp_V(itrial,icell) = mean(dataTrialStart(36:end,icell,v_ind(itrial)),1);
    end
end

baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);


cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsSelectZero_dir = intersect(dirSlctvCells,cellsPrefZero);
cellsSelectZero_ori = intersect(oriSlctvCells,cellsPrefZero);
cellsSelectZero = union(cellsSelectZero_dir,cellsSelectZero_ori);
cellsPrefRespZero = intersect(baselineStimRespIndex_V,cellsPrefZero);

cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
cellsSelectNinety_dir = intersect(dirSlctvCells,cellsPrefNinety);
cellsSelectNinety_ori = intersect(oriSlctvCells,cellsPrefNinety);
cellsSelectNinety = union(cellsSelectNinety_dir,cellsSelectNinety_ori);
cellsPrefRespNinety = intersect(baselineStimRespIndex_V,cellsPrefNinety);

nCells = size(cycDataDFoverF_cmlvNoTarget{7},2);
oriSlctvCellsAll = union(oriSlctvCells,dirSlctvCells);
notSlctvCells = setdiff([1:nCells],oriSlctvCellsAll);
notRespCells = setdiff([1:nCells],baselineStimRespIndex_V);

for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));   
    cycTrialOutcome{icyc} = trialOutcome(ind);
    cycDirectionDeg{icyc} = DirectionDeg(ind);
%     cycCatchDirectionDeg{icyc} = catchDirectionDeg(ind);
%     cycCatchTrialOutcome{icyc} = catchTrialOutcome(ind);
%     cycCatchCycle{icyc} = catchCycle(ind);
end
%%
dataavgtrials = squeeze(mean(dataDFoverF(:,:,V_cycInd),3));
drivencells = find(any(dataavgtrials >0.05,1));
%%
CYC = cycles(10);
cellgroup = 1:nCells;
downS = 1;    

tempdata = cycDataDFoverF_cmlvNoTarget{CYC};
Ltempdata = floor(size(tempdata,1)/downS)*downS;
tempdata = squeeze(mean(reshape(tempdata(1:Ltempdata,:,:), [downS floor(size(tempdata,1)/downS) size(tempdata,2) size(tempdata,3)]),1));
tempdata = permute(tempdata,[1,3,2]);

tempOutcome = cycTrialOutcome{CYC};
tempVind = intersect(cycV_ind{CYC},find(strcmp(tempOutcome,'success') == 1));
tempAVind = intersect(cycAV_ind{CYC},find(strcmp(tempOutcome,'success') == 1));

%%
cell1 = drivencells(5);
trialcorrV = corr(squeeze(tempdata(:,tempVind,cell1)));
trialcorrV = tril(trialcorrV,-1);
avgTrialcorrV = mean(mean(trialcorrV(trialcorrV ~= 0),2),1);

trialcorrAV = corr(squeeze(tempdata(:,tempAVind,cell1)));
trialcorrAV = tril(trialcorrAV,-1);
avgTrialcorrAV = mean(mean(trialcorrAV(trialcorrAV ~= 0),2),1);

figure;
colormap(brewermap([],'*RdBu'))
subplot(1,2,1);
imagesc(trialcorrV)
title([{['avg corr bw vis trials = ' num2str(avgTrialcorrV)]},{[num2str(length(tempVind)) ' trials']}])
colorbar 
caxis([-1 1]);
axis('square')
subplot(1,2,2);
imagesc(trialcorrAV)
title([{['avg corr bw aud trials = ' num2str(avgTrialcorrAV)]},{[num2str(length(tempAVind)) ' trials']}])
colorbar 
caxis([-1 1]);
axis('square')

%%
avgTrialcorrV = zeros(nCells,1);
avgTrialcorrAV = zeros(nCells,1);
for icell = 1:length(cellgroup)
    trialcorrV = corr(squeeze(tempdata(:,tempVind,icell)));
    trialcorrV = tril(trialcorrV,-1);
    avgTrialcorrV(icell,:) = mean(mean(trialcorrV(trialcorrV ~= 0),2),1);

    trialcorrAV = corr(squeeze(tempdata(:,tempAVind,icell)));
    trialcorrAV = tril(trialcorrAV,-1);
    avgTrialcorrAV(icell,:) = mean(mean(trialcorrAV(trialcorrAV ~= 0),2),1);
end

avgTrialcorrV_norm = avgTrialcorrV/max(avgTrialcorrV);
avgTrialcorrAV_norm = avgTrialcorrAV/max(avgTrialcorrV);

figure;
scatter(cellGroup,avgTrialcorrV_norm,'go')
hold on
scatter(cellgroup,avgTrialcorrAV_norm,'ko')