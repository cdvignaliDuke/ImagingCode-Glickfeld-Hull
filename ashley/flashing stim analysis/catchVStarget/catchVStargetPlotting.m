edit FlashingStim_dataSortedByCycle_combineDatasets.m

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
        dirTargetNorm{i} = [];
    elseif cTargetOn(ind(end),1)+40 > size(dataTC,1) 
        targetData = zeros(60,size(dataTC,2),length(ind));
        targetDataDF = zeros(60,size(dataTC,2),length(ind));
        targetDataDFoverF = zeros(60,size(dataTC,2),length(ind));
        targetBL = zeros(length(Dirs),size(dataTC,2),length(ind));
        targetNorm = zeros(60,size(dataTC,2),length(ind));
        for itrial = 1:length(ind)-1
            targetData(:,:,itrial) = dataTC(cTargetOn(ind(itrial))-19:cTargetOn(ind(itrial))+40,:);
            targetBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
            targetDataDF(:,:,itrial) = bsxfun(@minus, targetData(:,:,itrial), mean(targetBaseline,1));
            targetDataDFoverF(:,:,itrial) = bsxfun(@rdivide, targetDataDF(:,:,itrial), mean(targetBaseline,1));
            targetBL(i,:,itrial) = mean(targetDataDFoverF(17:21,:,itrial),1);
            targetNorm(:,:,itrial) = bsxfun(@minus,targetDataDFoverF(:,:,itrial),targetBL(i,:,itrial));
        end
        dirTargetData{i} = targetData;
        dirTargetDataDF{i} = targetDataDF;
        dirTargetDataDFoverF{i} = targetDataDFoverF;
        dirTargetInd{i} = ind;
        dirTargetNorm{i} = targetNorm;        
    else

        targetData = zeros(60,size(dataTC,2),length(ind));
        targetDataDF = zeros(60,size(dataTC,2),length(ind));
        targetDataDFoverF = zeros(60,size(dataTC,2),length(ind));
        targetBL = zeros(length(Dirs),size(dataTC,2),length(ind));
        targetNorm = zeros(60,size(dataTC,2),length(ind));
        for itrial = 1:length(ind)
            targetData(:,:,itrial) = dataTC(cTargetOn(ind(itrial))-19:cTargetOn(ind(itrial))+40,:);
            targetBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
            targetDataDF(:,:,itrial) = bsxfun(@minus, targetData(:,:,itrial), mean(targetBaseline,1));
            targetDataDFoverF(:,:,itrial) = bsxfun(@rdivide, targetDataDF(:,:,itrial), mean(targetBaseline,1));
            targetBL(i,:,itrial) = mean(targetDataDFoverF(17:21,:,itrial),1);
            targetNorm(:,:,itrial) = bsxfun(@minus,targetDataDFoverF(:,:,itrial),targetBL(i,:,itrial));
        end
        dirTargetData{i} = targetData;
        dirTargetDataDF{i} = targetDataDF;
        dirTargetDataDFoverF{i} = targetDataDFoverF;
        dirTargetInd{i} = ind;
        dirTargetNorm{i} = targetNorm;        
    end
    
end
%% sort trials by catch direction change

for i = 1:length(catchDirs)
%     ind = intersect(find(catchDirectionDeg == catchDirs(i)),find(strcmp(catchTrialOutcome,'CR') | strcmp(catchTrialOutcome,'FA')));
    ind = intersect(find(catchDirectionDeg == catchDirs(i)),find(strcmp(catchTrialOutcome,'FA')));
    if isempty(ind) == 1
            dirCatchData{i} = [];
            dirCatchDataDF{i} = [];
            dirCatchDataDFoverF{i} = [];
            dirCatchInd{i} = [];
            dirCatchNorm{i} = [];
    else
    catchData = zeros(60,size(dataTC,2),length(ind));
    catchDataDF = zeros(60,size(dataTC,2),length(ind));
    catchDataDFoverF = zeros(60,size(dataTC,2),length(ind));
    catchBL = zeros(length(catchDirs),size(dataTC,2),length(ind));
    catchNorm = zeros(60,size(dataTC,2),length(ind));
    
    for itrial = 1:length(ind)
        catchData(:,:,itrial) = dataTC((cCatchOn(ind(itrial)) - 19:cCatchOn(ind(itrial)) + 40),:);
        catchBaseline = dataTC(cLeverDown(ind(itrial))-29:cLeverDown(ind(itrial)),:);
        catchDataDF(:,:,itrial) = bsxfun(@minus,catchData(:,:,itrial),mean(catchBaseline,1));
        catchDataDFoverF(:,:,itrial) = bsxfun(@rdivide,catchDataDF(:,:,itrial),mean(catchBaseline,1));
        catchBL(i,:,itrial) = mean(catchDataDFoverF(17:21,:,itrial),1);
        catchNorm(:,:,itrial) = bsxfun(@minus,catchDataDFoverF(:,:,itrial),catchBL(i,:,itrial));
    end
    
    dirCatchData{i} = catchData;
    dirCatchDataDF{i} = catchDataDF;
    dirCatchDataDFoverF{i} = catchDataDFoverF;
    dirCatchInd{i} = ind;
    dirCatchNorm{i} = catchNorm;
    end
end

%% find mean response for each target and catch direction
meanTargetResp = zeros(60,size(dataTC,2),length(Dirs));
meanCatchResp = zeros(60,size(dataTC,2),length(catchDirs));
errTargetResp = zeros(60,size(dataTC,2),length(Dirs));
errCatchResp = zeros(60,size(dataTC,2),length(catchDirs));
for i = 1:length(Dirs)
    if i <= length(catchDirs)
        tempdataTarget = dirTargetNorm{i};
        tempdataCatch = dirCatchNorm{i};
        if ~isempty(tempdataTarget) 
            meanTargetResp(:,:,i) = mean(tempdataTarget,3);
            errTargetResp(:,:,i) = (std(tempdataTarget,0,3))./(sqrt(size(tempdataTarget,3)));
        end
        if ~isempty(tempdataCatch)
            meanCatchResp(:,:,i) = mean(tempdataCatch,3);
            errCatchResp(:,:,i) = (std(tempdataCatch,0,3))./(sqrt(size(tempdataCatch,3)));
        end
    else
        if ~isempty(tempdataTarget)
        tempdataTarget = dirTargetDataDFoverF{i};
        meanTargetResp(:,:,i) = mean(tempdataTarget,3);
        else
            meanTargetResp(:,:,i) = meanTargetResp(:,:,i);
        end
        
    end
end
        
%% plot
for i = 1: length(Dirs)
    legendInfoTarget{i} = [num2str(Dirs(i)) ' degrees'];
end
for i = 1: length(catchDirs)
    legendInfoCatch{i} = [num2str(catchDirs(i)) ' degrees'];
end

figure;
plot(squeeze(mean(meanTargetResp,2)))
vline(20,'c')
legend(legendInfoTarget)

figure;
plot(squeeze(mean(meanCatchResp,2)))
vline(20,'b')
legend(legendInfoCatch)

%% plot normalized catch resp
figure;
colors = brewermap(length(catchDirs)+2,'*Blues');
colors = colors(1:end-2,:);
% colors = brewermap(length(catchDirs)+15,'*RdPu');
% colorind = [3:2:length(catchDirs)+12];
% colors = colors(colorind,:);

for i = 1:length(catchDirs)
    tempdata = dirCatchNorm{i};
    if isempty(tempdata)
        tempdata = zeros(60,1);
        plot(nanmean(nanmean(tempdata,3),2),'Color',colors(i,:),'LineWidth',3)
        hold on
    else
        plot(nanmean(nanmean(tempdata,3),2),'Color',colors(i,:),'LineWidth',3)
        hold on
    end
end
legend(legendInfoCatch,'Location','NorthWest')
hold on
vline(20,'k')
vline(20+tooFastTime,'k:')
vline(20+maxReactTime,'k:')
title([mouse ' - ' date '; catch resp'])
% set(gca,'Color','k')


%%
cells =cellsPrefRespNinety;
orichange = catchDirs(end);

colorsC = brewermap(length(catchDirs)+2,'*Blues');
colorsC = colorsC(1:end-2,:);

colorsT = brewermap(length(Dirs)+15,'*YlGn');
colorindT = [3:2:length(Dirs)+12];
colorsT = colorsT(colorindT(1:length(Dirs)),:);


figure;
for i = 1:length(Dirs)
plot(mean(meanTargetResp(:,cells,i),2),'Color',colorsT(i,:),'LineWidth',3)
hold on
end
for i = 1:length(catchDirs)
plot(mean(meanCatchResp(:,cells,i),2),'Color',colorsC(i,:),'LineWidth',3)
hold on
end
legendInfoCombo = [legendInfoTarget legendInfoCatch];
legend(legendInfoCombo)

% plot(squeeze(mean(meanTargetResp(:,cells,find(Dirs == orichange)),2)),'g')
% hold on
% plot(squeeze(mean(meanCatchResp(:,cells,find(catchDirs == orichange)),2)),'b')
% hold on
vline(20,'c')
% title(['orichange = ' num2str(orichange) ';pref 90 deg'])
% legend({'target' 'catch'})

%%
cells =1:nCells;
cellgroupname = 'all cells';
orichange = catchDirs(end);

colorsC = brewermap(length(catchDirs)+2,'*Blues');
colorsC = colorsC(1:end-2,:);

colorsT = brewermap(length(Dirs)+15,'*YlGn');
colorindT = [3:2:length(Dirs)+12];
colorsT = colorsT(colorindT(1:length(Dirs)),:);

figure;
ind = find(ismember(Dirs,intersect(Dirs,catchDirs)));
for i = 1:length(ind)-1
    lineprops.col = {colorsT(1,:)};
    subplot(1,2,i)    
    mseb([],mean(meanTargetResp(:,cells,ind(i+1)),2),(std(meanTargetResp(:,cells,ind(i+1)),0,2)/sqrt(length(cells)))',lineprops,1);
    ylim([-0.01 0.03])
    legend({'target' 'catch'})
    hold on
end


for i = 1:length(catchDirs)-1
    lineprops.col = {colorsC(1,:)};
    subplot(1,2,i)
    mseb([],mean(meanCatchResp(:,cells,i+1),2),(std(meanCatchResp(:,cells,i+1),0,2)/sqrt(length(cells)))',lineprops,1)
    ylim([-0.01 0.03])
    legend({'target' 'catch'})
    title([num2str(catchDirs(i+1)) ' degrees'])
    hold on
end

suptitle([mouse '; ' date '; ' cellgroupname])
%%
cellGroup = 1:nCells;
% cellsSubgroup1 = find(ismember(cellGroup,intersect(cellGroup,cellsPrefZero)));
% cellsSubgroup2 = find(ismember(cellGroup,intersect(cellGroup,cellsPrefNinety)));
% cellsSubgroup3 = find(ismember(cellGroup,intersect(cellGroup,setdiff([1:nCells],cat(1,cellsPrefZero,cellsPrefNinety)))));
cellsSubgroup1 = cellsPrefZero;
cellsSubgroup2 = setdiff(cellGroup,cat(1,cellsPrefZero,cellsPrefNinety));
cellsSubgroup3 = cellsPrefNinety;
cellgroupname1 = 'pref 0 cells';
cellgroupname2 = 'other pref cells';
cellgroupname3 = 'pref90 cells';

% mean response to catch
cellCatchDirMean = zeros(60,nCells,length(catchDirs));
cellCatchDirMean_point = zeros(nCells,length(catchDirs));
cellCatchdirMean_diff = zeros(nCells,length(catchDirs));
for i = 1:length(catchDirs)
    tempdata = dirCatchNorm{i};
    if isempty(tempdata)
        cellCatchDirMean(:,:,i) = zeros(60,nCells);
        cellCatchDirMean_point(:,i) = squeeze(nanmean(cellCatchDirMean(25:29,:,i),1));
        cellCatchDirMean_diff(:,i) = squeeze(cellCatchDirMean(25,:,i)-cellCatchDirMean(29,:,i));
    else
        cellCatchDirMean(:,:,i) = nanmean(tempdata,3);
        cellCatchDirMean_point(:,i) = squeeze(nanmean(cellCatchDirMean(25:29,:,i),1));
        cellCatchDirMean_diff(:,i) = squeeze(cellCatchDirMean(25,:,i)-cellCatchDirMean(29,:,i));
    end
end
    
% mean response to target
cellTargetDirMean = zeros(60,nCells,length(Dirs));
cellTargetDirMean_point = zeros(nCells,length(Dirs));
cellTargetdirMean_diff = zeros(nCells,length(Dirs));
for i = 1:length(Dirs)
    if isempty(dirTargetNorm{i}) == 1
        cellTargetDirMean(:,:,i) = zeros(60,nCells);
        cellTargetDirMean_point(:,i) = zeros(nCells,1);
        cellTargetDirMean_diff(:,i) = zeros(nCells,1);
    else
        cellTargetDirMean(:,:,i) = nanmean(dirTargetNorm{i},3);
        cellTargetDirMean_point(:,i) = squeeze(nanmean(cellTargetDirMean(25:29,:,i),1));
        cellTargetDirMean_diff(:,i) = squeeze(cellTargetDirMean(25,:,i)-cellTargetDirMean(29,:,i));
    end
end
%% plot cell response to target by tuning preference (mean resp to target) (PLOT1)
figure;
colors = brewermap(length(catchDirs)+2,'*Blues');
colors = colors(1:end-2,:);
% colors = brewermap(length(Dirs)+15,'*RdPu');
% colorind = [3:2:length(Dirs)+12];
% colors = colors(colorind,:);

for idir = 1:length(catchDirs)
    subplot(3,3,idir)
    if idir == 1
        plot([1:length(cellsSubgroup1)],cellCatchDirMean_point(cellsSubgroup1,idir),'ko')
        title(cellgroupname1)
    else
        plot([1:length(cellsSubgroup1)],cellCatchDirMean_point(cellsSubgroup1,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
        title([num2str(catchDirs(idir)) ' deg'])
    end
%     ylim([-0.03 0.03])
    hold on
end

hold on


for idir = 1:length(catchDirs)
    subplot(3,3,idir+3)
    if idir == 1
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellCatchDirMean_point(cellsSubgroup2,idir),'ko')
        title(cellgroupname2)
    else
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellCatchDirMean_point(cellsSubgroup2,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
        
    end
%     ylim([-0.03 0.03])
    hold on
end

hold on


for idir = 1:length(catchDirs)
    subplot(3,3,idir+6)
    if idir == 1
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellCatchDirMean_point(cellsSubgroup3,idir),'ko')
        title(cellgroupname3)
    else
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellCatchDirMean_point(cellsSubgroup3,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
%     ylim([-0.03 0.03])
    hold on
end

% %% plot target response point for catch direction degrees (PLOT2)
% celltargetrespfig = figure; %Comment out if want overlay of catch and target responses
% colors = brewermap(length(Dirs),'*YlGn');
colors = brewermap(length(Dirs)+15,'*YlGn');
colorind = [1:2:length(Dirs)+10];
colors = colors(colorind,:);
ind = find(ismember(Dirs,catchDirs));

for idir = 1:length(catchDirs)
    subplot(3,3,idir)
    if idir == 1
        plot([1:length(cellsSubgroup1)],cellTargetDirMean_point(cellsSubgroup1,ind(idir)),'ko','MarkerFaceColor', 'k')
        title(cellgroupname1)
    else
        plot([1:length(cellsSubgroup1)],cellTargetDirMean_point(cellsSubgroup1,ind(idir)),'Color',colors(ind(idir),:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(ind(idir),:))
        title([num2str(catchDirs(idir)) ' deg'])
    end
%     ylim([-0.03 0.03])
    hold on
end

hold on

for idir = 1:length(catchDirs)
    subplot(3,3,idir+3)
    if idir == 1
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellTargetDirMean_point(cellsSubgroup2,ind(idir)),'ko','MarkerFaceColor', 'k')
        title(cellgroupname2)
    else
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellTargetDirMean_point(cellsSubgroup2,ind(idir)),'Color',colors(ind(idir),:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(ind(idir),:))
    end
%     ylim([-0.03 0.03])
    hold on
end

hold on

for idir = 1:length(catchDirs)
    subplot(3,3,idir+6)
    if idir == 1
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellTargetDirMean_point(cellsSubgroup3,ind(idir)),'ko','MarkerFaceColor', 'k')
        title(cellgroupname3)
    else
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellTargetDirMean_point(cellsSubgroup3,ind(idir)),'Color',colors(ind(idir),:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(ind(idir),:))
    end
%     ylim([-0.03 0.03])
    hold on
end

suptitle([mouse '; ' date '; gr = tar, blu = cat, bla = aud'])

%% is the 'supression response seen for most cells real?
for i = 1:length(Dirs)
    try
        suprCellsCatch{i} = find(cellCatchDirMean_point(:,i) < 0);
        suprCellsTarget{i} = find(cellTargetDirMean_point(:,i) < 0);
    catch
        suprCellsTarget{i} = find(cellTargetDirMean_point(:,i) < 0);
    end
end
    
cells = suprCellsTarget{end};
cellgroupname = 'cells supressed by 90 deg catch';
orichange = catchDirs(end);

colorsC = brewermap(length(catchDirs)+2,'*Blues');
colorsC = colorsC(1:end-2,:);

colorsT = brewermap(length(Dirs)+15,'*YlGn');
colorindT = [3:2:length(Dirs)+12];
colorsT = colorsT(colorindT(1:length(Dirs)),:);

figure;
ind = find(ismember(Dirs,intersect(Dirs,catchDirs)));
for i = 1:length(ind)-1
    lineprops.col = {colorsT(1,:)};
    subplot(1,2,i)    
    mseb([],mean(meanTargetResp(:,cells,ind(i+1)),2),(std(meanTargetResp(:,cells,ind(i+1)),0,2)/sqrt(length(cells)))',lineprops,1);
    ylim([-0.01 0.03])
    legend({'target' 'catch'})
    hold on
end


for i = 1:length(catchDirs)-1
    lineprops.col = {colorsC(1,:)};
    subplot(1,2,i)
    mseb([],mean(meanCatchResp(:,cells,i+1),2),(std(meanCatchResp(:,cells,i+1),0,2)/sqrt(length(cells)))',lineprops,1)
    ylim([-0.01 0.03])
    legend({'target' 'catch'})
    title([num2str(catchDirs(i+1)) ' degrees'])
    hold on
end

suptitle([mouse '; ' date '; ' cellgroupname])
%% plot cell response to target by tuning preference (mean resp to target)

cellcatchrespfig = figure;
colors = brewermap(length(catchDirs)+2,'*Blues');
colors = colors(1:end-2,:);
% colors = brewermap(length(Dirs)+15,'*RdPu');
% colorind = [3:2:length(Dirs)+12];
% colors = colors(colorind,:);

for idir = 1:length(catchDirs)
    if idir == 1
        plot([1:length(cellsSubgroup1)],cellCatchDirMean_diff(cellsSubgroup1,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([1:length(cellsSubgroup1)],cellCatchDirMean_diff(cellsSubgroup1,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

hold on
vline(length(cellsSubgroup1),'k')

for idir = 1:length(catchDirs)
    if idir == 1
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellCatchDirMean_diff(cellsSubgroup2,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellCatchDirMean_diff(cellsSubgroup2,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

hold on
vline(length(cellsSubgroup1)+length(cellsSubgroup2),'k')

for idir = 1:length(catchDirs)
    if idir == 1
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellCatchDirMean_diff(cellsSubgroup3,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellCatchDirMean_diff(cellsSubgroup3,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

legend(legendInfoCatch)
title([mouse ' - ' date ';mean diff to catch'])



%% scatter plot cells' response to catch vs. target 
cells = 1:nCells;
orichange = catchDirs(end);

% targetData = cellTargetDirMean_point(:,find(Dirs == orichange));
% catchData = cellCatchDirMean_point(:,find(catchDirs == orichange));
targetData = cellTargetDirMean_diff(:,find(Dirs == orichange));
catchData = cellCatchDirMean_diff(:,find(catchDirs == orichange));

colors = brewermap(4,'*Greys');
colors = colors(1:3,:);

scatter(targetData(cellsSubgroup1),catchData(cellsSubgroup1),100,colors(1,:),'filled')
hold on
scatter(targetData(cellsSubgroup2),catchData(cellsSubgroup2),100,colors(2,:),'filled')
hold on
scatter(targetData(cellsSubgroup3),catchData(cellsSubgroup3),100,colors(3,:),'filled')
hold on
legend({'zero' 'other' 'nintey'})
refline(1,0);
xlabel('target')
ylabel('catch')
% title([mouse ' - ' date ';mean target resp point'])
title([mouse ' - ' date ';mean target resp diff'])
% xlim([-0.01 0.03])
% ylim([-0.01 0.03])
axis('square')


%%
% trials = dirCatchInd{end};
% catchCyc = catchCycle(trials);
% trCyc = tCyclesOn(trials);
% 
% plotData = zeros(150,size(dataTC,2),length(trials));
% plotDataDF = zeros(150,size(dataTC,2),length(trials));
% plotDataDFoverF = zeros(150,size(dataTC,2),length(trials));
% for i = 1:length(trials)
%     plotData(:,:,i) = dataTC(cLeverDown(trials(i))-19:cLeverDown(trials(i))+130,:);
%     plotDataDF(:,:,i) = bsxfun(@minus, plotData(:,:,i), mean(plotData(1:20,:,i),1));
%     plotDataDFoverF(:,:,i) = bsxfun(@rdivide, plotDataDF(:,:,i), mean(plotData(1:20,:,i),1)); 
% end
% 
% for icell = 1:(length(cellsPrefNinety)/2)
% figure;
% for i = 1:length(trials)
%     subplot(3,4,i)
%     plot(squeeze(plotDataDFoverF(:,icell,i)))
% 	vline(20,'k')
%     hold on
%     icycC = catchCyc(i);
%     icycT = trCyc(i);
%     for imark = 1:cycles(icycT)-1
%         L = ((imark+1)*cycTime);
%         vline(L,'k:');
%         hold on
%     end
%     hold on
%     vline(((icycC+1)*cycTime),'b');
%     hold on
%     vline(((icycT+1)*cycTime),'c');
%     title(['catch cyc = ' num2str(icycC) ';target deg = ' num2str(DirectionDeg(trials(i)))])
% end
% end
