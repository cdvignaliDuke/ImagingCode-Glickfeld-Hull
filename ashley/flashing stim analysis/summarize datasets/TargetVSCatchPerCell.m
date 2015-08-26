% run(FlashingStim_dataSortedByCycle_combineDatasets.m)
%%
% cells selective for 90 deg, 0 deg, driven by baseline stim, "driven
% cells", non-selective cells, target driven cells
%%
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' date '\TargetandCatchAvgTraces'];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' date]);
        mkdir('TargetandCatchAvgTraces')
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(date,'TargetandCatchAvgTraces')
    end
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
%% find sets of cells
DirFolder = '006';
run('cellSets.m')

%% sort trials by direction change
for i = 1:length(Dirs)
%     ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'success')));    
    ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'success') | strcmp(trialOutcome,'ignore')));    
%         ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'ignore')));    
    if isempty(ind) == 1;
        dirTargetData{i} = [];
        dirTargetDataDF{i} = [];
        dirTargetDataDFoverF{i} = [];
        dirTargetInd{i} = ind;
%         dirTargetNorm{i} = [];
    elseif cTargetOn(ind(end),1)+40 > size(dataTC,1) 
        targetData = zeros(60,size(dataTC,2),length(ind));
        targetDataDF = zeros(60,size(dataTC,2),length(ind));
        targetDataDFoverF = zeros(60,size(dataTC,2),length(ind));
%         targetBL = zeros(length(Dirs),size(dataTC,2),length(ind));
%         targetNorm = zeros(60,size(dataTC,2),length(ind));
        for itrial = 1:length(ind)-1
            targetData(:,:,itrial) = dataTC(cTargetOn(ind(itrial))-19:cTargetOn(ind(itrial))+40,:);
            targetBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
            targetDataDF(:,:,itrial) = bsxfun(@minus, targetData(:,:,itrial), mean(targetBaseline,1));
            targetDataDFoverF(:,:,itrial) = bsxfun(@rdivide, targetDataDF(:,:,itrial), mean(targetBaseline,1));
%             targetBL(i,:,itrial) = mean(targetDataDFoverF(17:21,:,itrial),1);
%             targetNorm(:,:,itrial) = bsxfun(@minus,targetDataDFoverF(:,:,itrial),targetBL(i,:,itrial));
        end
        dirTargetData{i} = targetData;
        dirTargetDataDF{i} = targetDataDF;
        dirTargetDataDFoverF{i} = targetDataDFoverF;
        dirTargetInd{i} = ind;
%         dirTargetNorm{i} = targetNorm;        
    else
        targetData = zeros(60,size(dataTC,2),length(ind));
        targetDataDF = zeros(60,size(dataTC,2),length(ind));
        targetDataDFoverF = zeros(60,size(dataTC,2),length(ind));
%         targetBL = zeros(length(Dirs),size(dataTC,2),length(ind));
%         targetNorm = zeros(60,size(dataTC,2),length(ind));
        for itrial = 1:length(ind)
            targetData(:,:,itrial) = dataTC(cTargetOn(ind(itrial))-19:cTargetOn(ind(itrial))+40,:);
            targetBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
            targetDataDF(:,:,itrial) = bsxfun(@minus, targetData(:,:,itrial), mean(targetBaseline,1));
            targetDataDFoverF(:,:,itrial) = bsxfun(@rdivide, targetDataDF(:,:,itrial), mean(targetBaseline,1));
%             targetBL(i,:,itrial) = mean(targetDataDFoverF(17:21,:,itrial),1);
%             targetNorm(:,:,itrial) = bsxfun(@minus,targetDataDFoverF(:,:,itrial),targetBL(i,:,itrial));
        end
        dirTargetData{i} = targetData;
        dirTargetDataDF{i} = targetDataDF;
        dirTargetDataDFoverF{i} = targetDataDFoverF;
        dirTargetInd{i} = ind;
%         dirTargetNorm{i} = targetNorm;        
    end
    
end

tempdata = nanmean(dirTargetDataDFoverF{end},3);
if isempty(tempdata) == 1
    tempdata = nanmean(dirTargetDataDFoverF{end-1},3);
end
meanBefore = nanmean(tempdata(1:23,:),1);
meanAfter = nanmean(tempdata(24:end,:),1);
meanSub = meanAfter-meanBefore;
targetdrivencells = find(meanSub >0.05);
%% sort trials by catch direction change

for i = 1:length(catchDirs)
    ind = intersect(find(catchDirectionDeg == catchDirs(i)),find(strcmp(catchTrialOutcome,'CR') | strcmp(catchTrialOutcome,'FA')));
%     ind = intersect(find(catchDirectionDeg == catchDirs(i)),find(strcmp(catchTrialOutcome,'CR')));
    if isempty(ind) == 1
            dirCatchData{i} = [];
            dirCatchDataDF{i} = [];
            dirCatchDataDFoverF{i} = [];
            dirCatchInd{i} = [];
%             dirCatchNorm{i} = [];
    else
    catchData = zeros(60,size(dataTC,2),length(ind));
    catchDataDF = zeros(60,size(dataTC,2),length(ind));
    catchDataDFoverF = zeros(60,size(dataTC,2),length(ind));
%     catchBL = zeros(length(catchDirs),size(dataTC,2),length(ind));
%     catchNorm = zeros(60,size(dataTC,2),length(ind));
    
    for itrial = 1:length(ind)
        catchData(:,:,itrial) = dataTC((cCatchOn(ind(itrial)) - 19:cCatchOn(ind(itrial)) + 40),:);
        catchBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
        catchDataDF(:,:,itrial) = bsxfun(@minus,catchData(:,:,itrial),mean(catchBaseline,1));
        catchDataDFoverF(:,:,itrial) = bsxfun(@rdivide,catchDataDF(:,:,itrial),mean(catchBaseline,1));
%         catchBL(i,:,itrial) = mean(catchDataDFoverF(17:21,:,itrial),1);
%         catchNorm(:,:,itrial) = bsxfun(@minus,catchDataDFoverF(:,:,itrial),catchBL(i,:,itrial));
    end
    
    dirCatchData{i} = catchData;
    dirCatchDataDF{i} = catchDataDF;
    dirCatchDataDFoverF{i} = catchDataDFoverF;
    dirCatchInd{i} = ind;
%     dirCatchNorm{i} = catchNorm;
    end
end

%% ***************************************************
cells =1:nCells;
cellgroupname = 'all cells';
figBaseName = 'avgDFoverFpretarget_allcells_all';

%***************************

%% find mean response for each target and catch direction
meanTargetRespNorm_cells = zeros(60,length(cells),length(Dirs));
meanTargetRespNorm = zeros(60,length(Dirs));
meanCatchRespNorm_cells = zeros(60,length(cells),length(catchDirs));
meanCatchRespNorm = zeros(60,length(catchDirs));
errTargetResp_cells = zeros(60,length(cells),length(Dirs));
errTargetResp = zeros(60,length(Dirs));
errCatchResp_cells = zeros(60,length(cells),length(catchDirs));
errCatchResp = zeros(60,length(catchDirs));
nTrialsTarget = zeros(1,length(Dirs));
nTrialsCatch = zeros(1,length(catchDirs));

for i = 1:length(Dirs)
    tempdataTarget = dirTargetDataDFoverF{i};
    if ~isempty(tempdataTarget)
        meanTargetRespNorm_cells(:,:,i) = mean(tempdataTarget(:,cells,:),3);
        meanTargetRespNorm(:,i) = mean(mean(tempdataTarget(:,cells,:),3),2);
        errTargetResp_cells(:,:,i) = (std(tempdataTarget(:,cells,:),0,3))./(sqrt(size(tempdataTarget(:,cells,:),3)));
        errTargetResp(:,i) = (squeeze(std(meanTargetRespNorm_cells(:,:,i),0,2)))./(sqrt(size(tempdataTarget(:,cells,:),2)));
        nTrialsTarget(:,i) = size(tempdataTarget,3);
    end
end
meanTargetRespNorm_cells = bsxfun(@minus,meanTargetRespNorm_cells,mean(meanTargetRespNorm_cells(17:21,:,:),1));
meanTargetRespNorm = bsxfun(@minus,meanTargetRespNorm,mean(meanTargetRespNorm(17:21,:),1));
errTargetResp_cells = bsxfun(@minus,errTargetResp_cells,mean(meanTargetRespNorm_cells(17:21,:,:),1));
errTargetResp = bsxfun(@minus,errTargetResp,mean(meanTargetRespNorm(17:21,:),1));

for i = 1:length(catchDirs)
    tempdataCatch = dirCatchDataDFoverF{i};
    if ~isempty(tempdataCatch)
        meanCatchRespNorm_cells(:,:,i) = mean(tempdataCatch(:,cells,:),3);
        meanCatchRespNorm(:,i) = mean(mean(tempdataCatch(:,cells,:),3),2);
        errCatchResp_cells(:,:,i) = (std(tempdataCatch(:,cells,:),0,3))./(sqrt(size(tempdataCatch(:,cells,:),3)));
        errCatchResp(:,i) = (squeeze(std(meanCatchRespNorm_cells(:,:,i),0,2)))./(sqrt(size(tempdataCatch(:,cells,:),2)));
        nTrialsCatch(:,i) = size(tempdataCatch,3);
    end
end
meanCatchRespNorm_cells = bsxfun(@minus,meanCatchRespNorm_cells,mean(meanCatchRespNorm_cells(17:21,:,:),1));
meanCatchRespNorm = bsxfun(@minus,meanCatchRespNorm,mean(meanCatchRespNorm(17:21,:),1));
errCatchResp_cells = bsxfun(@minus,errCatchResp_cells,mean(meanCatchRespNorm_cells(17:21,:,:),1));
errCatchResp = bsxfun(@minus,errCatchResp,mean(meanCatchRespNorm(17:21,:),1));

%%
meanTargetRespNorm_cellsDiff = squeeze(mean(meanTargetRespNorm_cells(23:27,:,:),1)) - squeeze(mean(meanTargetRespNorm_cells(18:22,:,:),1));
meanCatchRespNorm_cellsDiff = squeeze(mean(meanCatchRespNorm_cells(23:27,:,:),1)) - squeeze(mean(meanCatchRespNorm_cells(18:22,:,:),1));

figure;
scatter(meanTargetRespNorm_cellsDiff(:,end),meanCatchRespNorm_cellsDiff(:,end))
hold on
refline([1,0]);
xlabel('target')
ylabel('catch')