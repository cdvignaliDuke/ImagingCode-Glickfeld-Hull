iexp = 11;
%%
ialign = 1;
run('divideupdatabyalignment.m')
%% find sets of cells
DirFolder = expt(iexp).dirtuning;
run('cellSets.m')
%%
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' date '\PreTargetAvgTraces'];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' date]);
        mkdir('PreTargetAvgTraces')
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(date,'PreTargetAvgTraces')
    end
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

%% sort trials by direction change
for i = 1:length(Dirs)
    ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'success')));    
%     ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'success') | strcmp(trialOutcome,'ignore')));    
%         ind = intersect(find(DirectionDeg == Dirs(i)),find(strcmp(trialOutcome,'ignore')));    
    if isempty(ind) == 1;
        dirTargetData{i} = [];
        dirTargetDataDF{i} = [];
        dirTargetDataDFoverF{i} = [];
        dirTargetInd{i} = ind;
%         dirTargetNorm{i} = [];
    elseif cTargetOn(1,ind(end))+40 > size(dataTC,1) 
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
%     ind = intersect(find(catchDirectionDeg == catchDirs(i)),find(strcmp(catchTrialOutcome,'CR') | strcmp(catchTrialOutcome,'FA')));
    ind = intersect(find(catchDirectionDeg == catchDirs(i)),find(strcmp(catchTrialOutcome,'FA')));
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
cells = 1:nCells;
cellgroupname = 'all cells';
figBaseName = 'avgDFoverFpretarget_allcells_successignore';

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

%% legend and color coding
for i = 1: length(catchDirs)
    legendInfoCatch{i} = [num2str(catchDirs(i)) ' degrees; ' num2str(nTrialsCatch(i)) ' trials'];
end
for i = 1: length(Dirs)
    legendInfoTarget{i} = [num2str(Dirs(i)) ' degrees; ' num2str(nTrialsTarget(i)) ' trials'];
end
legendInfoTarget{1} = 'aud., 0 degrees';

colorsC = brewermap(length(catchDirs)+2,'*Blues');
colorsC = colorsC(1:end-2,:);

colorsT = brewermap(length(Dirs)+15,'*YlGn');
colorindT = [3:2:length(Dirs)+12];
colorsT = colorsT(colorindT(1:length(Dirs)),:);
%%
meanTargetRespNorm1 = meanTargetRespNorm;
errTargetResp1 = errTargetResp;
cellgroupname1 = cellgroupname;
% meanTargetRespNorm2 = meanTargetRespNorm;
% errTargetResp2 = errTargetResp;
% cellgroupname2 = cellgroupname;
% meanTargetRespNorm3 = meanTargetRespNorm;
% errTargetResp3 = errTargetResp;
% cellgroupname3 = cellgroupname;
% meanTargetRespNorm4 = meanTargetRespNorm;
% errTargetResp4 = errTargetResp;
% cellgroupname4 = cellgroupname;
%% plot all target responses
% figName = '_avgresp2target';
% figure;
% for i = 1:length(Dirs)
%     if i == 1
%         lineprops.col = {'k'};
%     else
%         lineprops.col = {colorsT(i,:)};
%     end
% %     subplot(1,2,i)    
%     mseb([-19:40],meanTargetRespNorm(:,i),errTargetResp(:,i)',lineprops,1);
% %     ylim([-0.01 0.03])    
%     hold on
% end
% legend(legendInfoTarget,'Location', 'NorthWest')
% title({[mouse '; ' date]; 'target resp'; cellgroupname})
% print([fnout ['\' figBaseName figName '.pdf']], '-dpdf')
%% plot all target responses
figName = '_avgresp2target';
figure;
for i = 1:length(Dirs)
    if i == 1
        lineprops.col = {'k'};
    else
        lineprops.col = {colorsT(i,:)};
    end
%     subplot(1,2,i)    
    mseb([-19:40],meanTargetRespNorm(:,i),errTargetResp(:,i)',lineprops,1);
%     ylim([-0.01 0.03])    
    hold on
end
legend(legendInfoTarget,'Location', 'NorthWest')
title({[mouse '; ' date]; 'target resp'; cellgroupname})
print([fnout ['\' figBaseName figName '.pdf']], '-dpdf')


%% 
figure;
subplot(2,2,1)
for i = 1:length(Dirs)
    if i == 1
        lineprops.col = {'k'};
    else
        lineprops.col = {colorsT(i,:)};
    end
%     subplot(1,2,i)    
    mseb([-19:40],meanTargetRespNorm1(:,i),errTargetResp1(:,i)',lineprops,1);
%     ylim([-0.01 0.03])    
    hold on
end
legend(legendInfoTarget,'Location', 'NorthWest')
title({[mouse '; ' date]; 'target resp'; cellgroupname1})

subplot(2,2,2)
for i = 1:length(Dirs)
    if i == 1
        lineprops.col = {'k'};
    else
        lineprops.col = {colorsT(i,:)};
    end
%     subplot(1,2,i)    
    mseb([-19:40],meanTargetRespNorm2(:,i),errTargetResp2(:,i)',lineprops,1);
%     ylim([-0.01 0.03])    
    hold on
end
title({[mouse '; ' date]; 'target resp'; cellgroupname2})

subplot(2,2,3)
for i = 1:length(Dirs)
    if i == 1
        lineprops.col = {'k'};
    else
        lineprops.col = {colorsT(i,:)};
    end
%     subplot(1,2,i)    
    mseb([-19:40],meanTargetRespNorm3(:,i),errTargetResp3(:,i)',lineprops,1);
%     ylim([-0.01 0.03])    
    hold on
end
title({[mouse '; ' date]; 'target resp'; cellgroupname3})

subplot(2,2,4)
for i = 1:length(Dirs)
    if i == 1
        lineprops.col = {'k'};
    else
        lineprops.col = {colorsT(i,:)};
    end
%     subplot(1,2,i)    
    mseb([-19:40],meanTargetRespNorm4(:,i),errTargetResp4(:,i)',lineprops,1);
%     ylim([-0.01 0.03])    
    hold on
end
title({[mouse '; ' date]; 'target resp'; cellgroupname4})

print([fnout ['\' 'combinetargetresp' '.pdf']], '-dpdf')
%% plot all catch responses
figure;
% subplot(2,2,1)
for i = 1:length(catchDirs)
    if i == 1
        lineprops.col = {'k'};
    else
        lineprops.col = {colorsC(i,:)};
    end
%     subplot(2,2,i)
    mseb([-19:40],meanCatchRespNorm(:,i),errCatchResp(:,i)',lineprops,1);
%     ylim([-0.01 0.03])
    hold on
end
legend(legendInfoCatch,'Location', 'NorthWest')
title({[mouse '; ' date]; 'catch resp'; cellgroupname1})
%% overlay catch and target responses
figName = 'FACR_catchVStargetresp';
figure;
ind = find(ismember(Dirs,intersect(Dirs,catchDirs)));
for i = 1:length(ind)-1
    lineprops.col = {colorsT(i,:)};
    subplot(1,2,i)
    mseb([-19:40],meanTargetRespNorm(:,ind(i+1)),errTargetResp(:,ind(i+1))',lineprops,1);
%     ylim([-0.01 0.05])
    legend({'target' 'catch'})
    hold on
end


for i = 1:length(catchDirs)-1
    lineprops.col = {colorsC(1,:)};
    subplot(1,2,i)
    mseb([-19:40],meanCatchRespNorm(:,i+1),errCatchResp(:,i+1)',lineprops,1);
%     ylim([-0.01 0.05])
    legend({'target' 'catch'})
    title([num2str(catchDirs(i+1)) ' degrees'])
    hold on
end

suptitle([mouse '; ' date '; ' cellgroupname])
print([fnout ['\' figBaseName figName '.pdf']], '-dpdf')
