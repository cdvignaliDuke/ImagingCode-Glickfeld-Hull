edit FlashingStim_dataSortedByCycle_combineDatasets.m

for icyc = 1:length(cycles)
    ind = find(tCyclesOn == cycles(icyc));   
    cycTrialOutcome{icyc} = trialOutcome(ind);
    cycDirectionDeg{icyc} = DirectionDeg(ind);
end

% Success_ind = find(strcmp('success',trialOutcome));
Dirs = unique(DirectionDeg);
% DirectionDeg_success_ind = DirectionDeg(Success_ind);
% 
% load('cellSelectivityIndices.mat')
% cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
% cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
% cellsSlctvZero = intersect(cellsPrefZero,dirSlctvCells);

%% color-code cells by orientation preference
DirFolder = '006';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, DirFolder);
cd(fileSave);
% load('oriTuningPreferences.mat')
load('TuningPreferences.mat')

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

%% sort cell trials by target grating orientation (average each) (successes only)
for icyc = 1:length(cycles)
    tempdata = cycDataDFoverF{icyc};
    cycdirs = cycDirectionDeg{icyc};
    cycoutcome = cycTrialOutcome{icyc};
    cycsuccess = find(strcmp(cycoutcome,'success'));
    cycdirmean = zeros(size(tempdata,1),size(tempdata,2),length(Dirs));
    for i = 1: length(Dirs)
        ind = intersect(cycsuccess, find(cycdirs == Dirs(i)));
        if isempty(ind)
            cycdirmean = [];
        else
            cycdirmean(:,:,i) = mean(tempdata(:,:,ind),3);
        end
    end
    cycTargetDirMean{icyc} = cycdirmean;
end

%% plot target response, ea cycle
figure;
for icyc = 1:length(cycles)
    tempdata = cycTargetDirMean{icyc};
    tempdatamean = squeeze(mean(tempdata,2));
    
    subplot(2,5,icyc)
    plot(tempdatamean(:,end))
    hold on
    vline(30,'k');
    hold on
    for i = 1:cycles(icyc)-1
        L = (i*cycTime)+30;
        vline(L,'k:');
        hold on
    end
    vline((cycles(icyc)*cycTime+30),'c');
    hold on
end
    
hold on
for i = 1: length(Dirs)
    legendInfo{i} = [num2str(Dirs(i)) ' degrees'];
end
legend(legendInfo,'Location','SouthEast')

%% Normalize target response
for icyc = 1:length(cycles)
    tempdata = cycTargetDirMean{icyc};
    targetindex = cycles(icyc)*cycTime+30;
    if isempty(tempdata) == 1
        targetBaselineF(icyc,:,:) = NaN(1,nCells,length(Dirs));
    else
        targetBaselineF(icyc,:,:) = mean(tempdata(targetindex-3:targetindex+2,:,:),1);
    end
end

for icyc = 1:length(cycles)
    tempdata = cycTargetDirMean{icyc};
    targetindex = cycles(icyc)*cycTime+30;
    if isempty(tempdata) == 1
        dirtraces(:,:,:,icyc) = NaN(40,nCells,length(Dirs));
    else
        tempdata = tempdata(targetindex-10:targetindex+29,:,:);
        tempdatasub = bsxfun(@minus,tempdata,targetBaselineF(icyc,:,:));
        for idir = 1:length(Dirs)
            dirtraces(:,:,idir,icyc) = tempdatasub(:,:,idir);
        end
    end
end

%%
targetDirMean = squeeze(nanmean(nanmean(dirtraces,4),2));

figure;
plot(targetDirMean)
hold on
legend(legendInfo,'Location','NorthWest')
hold on
vline(10,'k')
        

%% 
% set cell types, length of trial

cellGroup = 1:nCells;
% cellsSubgroup1 = find(ismember(cellGroup,intersect(cellGroup,cellsPrefZero)));
% cellsSubgroup2 = find(ismember(cellGroup,intersect(cellGroup,cellsPrefNinety)));
% cellsSubgroup3 = find(ismember(cellGroup,intersect(cellGroup,setdiff([1:nCells],cat(1,cellsPrefZero,cellsPrefNinety)))));
cellsSubgroup1 = cellsPrefZero;
cellsSubgroup2 = setdiff(cellGroup,cat(1,cellsPrefZero,cellsPrefNinety));
cellsSubgroup3 = cellsPrefNinety;

% mean response to target
cellTargetDirMean = nanmean(dirtraces,4);
cellTargetDirMean_point = squeeze(nanmean(cellTargetDirMean(21:25,:,:),1));

%% plot cell response to target by tuning preference

figure;
colors = brewermap(length(Dirs),'*Spectral');
% colors = brewermap(length(Dirs)+15,'*Greens');
% colorind = [3:2:length(Dirs)+12];
% colors = colors(colorind,:);

for idir = 1:length(Dirs)
    if idir == 1
        plot([1:length(cellsSubgroup1)],cellTargetDirMean_point(cellsSubgroup1,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([1:length(cellsSubgroup1)],cellTargetDirMean_point(cellsSubgroup1,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

hold on
vline(length(cellsSubgroup1),'k')

for idir = 1:length(Dirs)
    if idir == 1
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellTargetDirMean_point(cellsSubgroup2,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([length(cellsSubgroup1)+1:length(cellsSubgroup1)+length(cellsSubgroup2)],cellTargetDirMean_point(cellsSubgroup2,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

hold on
vline(length(cellsSubgroup1)+length(cellsSubgroup2),'k')

for idir = 1:length(Dirs)
    if idir == 1
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellTargetDirMean_point(cellsSubgroup3,idir),'ko','MarkerFaceColor', 'k')
    else
        plot([length(cellsSubgroup1)+length(cellsSubgroup2)+1:length(cellsSubgroup1)+length(cellsSubgroup2)+length(cellsSubgroup3)],cellTargetDirMean_point(cellsSubgroup3,idir),'Color',colors(idir,:),'LineStyle','none','Marker','o','MarkerFaceColor', colors(idir,:))
    end
    hold on
end

legend(legendInfo)
title([mouse ' - ' date])


%%
% L = 60;
% data_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind));
% dF_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind)); 
% dFoverF_aroundTarget = zeros(L,size(dataTC,2),length(Success_ind));
% for i = 1:length(Success_ind)
%     trial = Success_ind(i);
%     data_aroundTarget(:,:,i) = dataTC(cTargetOn(trial)-(L/2):cTargetOn(trial)+((L/2)-1),:);
%     dF_aroundTarget(:,:,i) = bsxfun(@minus, data_aroundTarget(:,:,i), mean(data_aroundTarget((L/2)-5:(L/2),:,i),1));
%     dFoverF_aroundTarget(:,:,i) = bsxfun(@rdivide, dF_aroundTarget(:,:,i),mean(data_aroundTarget((L/2)-5:(L/2),:,i),1));
% end
% 
% 
% V_Success_ind = find(ismember(Success_ind,intersect(Success_ind,V_ind)));
% 
% V_targetDirAvg = zeros(L,length(Dirs));
% %% cells pref 90
% for i = 1:length(Dirs)
%     ind = find(DirectionDeg_success_ind == Dirs(i));
%     ind = intersect(V_Success_ind,ind);
%     V_targetDirAvg90(:,i) = mean(mean(dFoverF_aroundTarget(:,cellsPrefNinety,ind),3),2);
% end
% 
% for i = 1:length(Dirs)
%     direction = Dirs(i);
%     legendinfo{i} = num2str(direction);
% end
% 
% cMap = colormap(winter(length(Dirs)));
% figure;
% for i = 1:length(Dirs)
%     plot(V_targetDirAvg90((20:end),i),'Color',cMap(i,:))
%     hold on
% end
% vline(10,'c')
% hold on
% legend(legendinfo)
% ylabel('dF/F')
% xlabel('frames')
% title('average target response, cells pref 90')
% 
% %% cells pref 0
% for i = 1:length(Dirs)
%     ind = find(DirectionDeg_success_ind == Dirs(i));
%     ind = intersect(V_Success_ind,ind);
%     V_targetDirAvg0(:,i) = mean(mean(dFoverF_aroundTarget(:,cellsPrefZero,ind),3),2);
% end
% 
% for i = 1:length(Dirs)
%     direction = Dirs(i);
%     legendinfo{i} = num2str(direction);
% end
% 
% cMap = colormap(winter(length(Dirs)));
% figure;
% for i = 1:length(Dirs)
%     plot(V_targetDirAvg0((20:end),i),'Color',cMap(i,:))
%     hold on
% end
% vline(10,'c')
% hold on
% legend(legendinfo)
% ylabel('dF/F')
% xlabel('frames')
% title('average target response, cells pref 0')
%     
% 
% 
