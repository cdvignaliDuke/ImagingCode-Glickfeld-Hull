% run(FlashingStim_dataSortedByCycle_combineDatasets.m)
%%
% cells selective for 90 deg, 0 deg, driven by baseline stim, "driven
% cells", non-selective cells
%%
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' date '\PreTargetRespPerCell'];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' date]);
        mkdir('PreTargetRespPerCell')
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(date,'PreTargetRespPerCell')
    end
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
%% find sets of cells
DirFolder = '006';
run('cellSets.m')
%%
cells = oriSlctvCellsAll;
cellgroupname = 'ori or dir selective';
figName = 'scatteravgDFoverFpretarget_orislctv_success';

%%
cellsSubgroup1 = find(ismember(cells,intersect(cells,cellsPrefZero)));
cellsSubgroup2 = find(ismember(cells,intersect(cells,cellsPrefNinety)));
cellsSubgroup3 = find(ismember(cells,intersect(cells,setdiff(cells,cat(1,cellsPrefZero,cellsPrefNinety)))));

%%
% scatter plot avg end of trial dF/F aud vs. vis trials
figure;
start = 1;
for icyc = 4:length(cycles)
    tempdata = cycDataDFoverF_cmlvNoTarget{icyc};
    

    tempdata = tempdata(end-29:end,:,:); %%%%% CHANGE THIS TO SELECT FRAMES TO ANALYZE
    
    V_ind = cycV_ind{icyc};
    AV_ind = cycAV_ind{icyc};
    V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    V_avg = squeeze(mean(mean(tempdata(:,cells,V_ind),3),1));
    AV_avg = squeeze(mean(mean(tempdata(:,cells,AV_ind),3),1));
    errbar_V = (std(squeeze(mean(tempdata(:,cells,V_ind),1)),[],2))./(sqrt(size(tempdata(:,cells,V_ind),3)));
    errbar_AV = (std(squeeze(mean(tempdata(:,cells,AV_ind),1)),[],2))./(sqrt(size(tempdata(:,cells,AV_ind),3)));
    subplot(2,4,start)
    ploterr(V_avg(:,cellsSubgroup1),AV_avg(:,cellsSubgroup1),errbar_V(cellsSubgroup1,:),errbar_AV(cellsSubgroup1,:),'go');
    hold on
    ploterr(V_avg(:,cellsSubgroup2),AV_avg(:,cellsSubgroup2),errbar_V(cellsSubgroup2,:),errbar_AV(cellsSubgroup2,:),'co');
    hold on
    ploterr(V_avg(:,cellsSubgroup3),AV_avg(:,cellsSubgroup3),errbar_V(cellsSubgroup3,:),errbar_AV(cellsSubgroup3,:),'ko');
    hold on
    scatter(mean(V_avg(:,cellsSubgroup1)),mean(AV_avg(:,cellsSubgroup1)),'g','filled');
    hold on
    scatter(mean(V_avg(:,cellsSubgroup2)),mean(AV_avg(:,cellsSubgroup2)),'c','filled');
    hold on
    scatter(mean(V_avg(:,cellsSubgroup3)),mean(AV_avg(:,cellsSubgroup3)),'k','filled');
	hold on
    refline(1,0);
    xlabel('visual')
    ylabel('auditory')
    title({[num2str(length(V_ind)) 'visual & ']; [num2str(length(AV_ind)) 'aud trials; ' num2str(icyc) ' cyc']});
    xlim([-0.1 0.3])
    ylim([-0.1 0.3])
    axis('square')
    start = start+1;
end

suptitle([mouse '; ' date '; ' cellgroupname])

print([fnout ['\' figName '.pdf']], '-dpdf')


%% ****WIP below here****
%% ratio of first to last base stim resp
figure;
for icyc = 1:length(cycles)
tempdata = cycDataDFoverF_cmlvNoTarget{icyc};
V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));

tempdataFirst = squeeze(mean(tempdata(33:37,cells,:),1)) - squeeze(mean(tempdata(28:32,cells,:),1));
lastCycInd = 30+ (cycTime*(icyc-1));
tempdataLast = squeeze(mean(tempdata(lastCycInd+3:lastCycInd+7,cells,:),1)) - squeeze(mean(tempdata(lastCycInd-2:lastCycInd+2,cells,:),1));

subplot(3,4,icyc)
scatter(mean(tempdataLast(cellsSubgroup1,V_cycInd),2),mean(tempdataFirst(cellsSubgroup1,V_cycInd),2),'g')
hold on
scatter(mean(tempdataLast(cellsSubgroup2,V_cycInd),2),mean(tempdataFirst(cellsSubgroup2,V_cycInd),2),'c')
hold on
scatter(mean(tempdataLast(cellsSubgroup3,V_cycInd),2),mean(tempdataFirst(cellsSubgroup3,V_cycInd),2),'k')
hold on
refline(1,0);
xlabel('Last')
ylabel('First')
title([num2str(length(V_cycInd)) 'visual trials; ' num2str(icyc) ' cyc']);
xlim([-0.02 0.03])
ylim([-0.02 0.03])
axis('square')
end

figure;
for icyc = 1:length(cycles)
tempdata = cycDataDFoverF_cmlvNoTarget{icyc};
V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));

tempdataFirst = squeeze(mean(tempdata(33:37,cells,:),1)) - squeeze(mean(tempdata(28:32,cells,:),1));
lastCycInd = 30+ (cycTime*(icyc-1));
tempdataLast = squeeze(mean(tempdata(lastCycInd+3:lastCycInd+7,cells,:),1)) - squeeze(mean(tempdata(lastCycInd-2:lastCycInd+2,cells,:),1));

% maxFirst = max(max(tempdataFirst,[],2),[],1);
% tempdataFirstNorm = bsxfun(@rdivide,tempdataFirst, maxFirst);
% tempdataLastNorm = bsxfun(@rdivide,tempdataLast, maxFirst);

tempdataRatio = tempdataFirst./tempdataLast;
maxRatio = max(max(tempdataRatio,[],2),[],1);
tempdataRatioN = bsxfun(@rdivide,tempdataRatio, maxRatio);

subplot(3,4,icyc)
scatter(mean(tempdataRatioN(cellsSubgroup1,V_cycInd),2),mean(tempdataRatioN(cellsSubgroup1,AV_cycInd),2),'g')
hold on
scatter(mean(tempdataRatioN(cellsSubgroup2,V_cycInd),2),mean(tempdataRatioN(cellsSubgroup2,AV_cycInd),2),'c')
hold on
scatter(mean(tempdataRatioN(cellsSubgroup3,V_cycInd),2),mean(tempdataRatioN(cellsSubgroup3,AV_cycInd),2),'k')
hold on
refline(1,0);
xlabel('Visual')
ylabel('Auditory')
title([num2str(length(V_cycInd)) 'visual trials; ' num2str(icyc) ' cyc']);
% xlim([-0.02 0.03])
% ylim([-0.02 0.03])
axis('square')
end