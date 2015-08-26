% run(FlashingStim_dataSortedByCycle_combineDatasets.m)
%%
% oriselectivecells, driven by baseline stim, "driven
% cells", non-selective cells, target driven cells
%%
fnout = ['Z:\Analysis\' mouse '\two-photon imaging\' date '\PreTargetRespReliability'];
try
    cd(fnout)
catch
    try
        cd(['Z:\Analysis\' mouse '\two-photon imaging\' date]);
        mkdir('PreTargetRespReliability')
    catch
        cd(['Z:\Analysis\' mouse '\two-photon imaging\']);
        mkdir(date,'PreTargetRespReliability')
    end
end

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
%% find sets of cells
DirFolder = '006';
run('cellSets.m')
%%
cells = drivencells;
cellgroupname = 'drivencells';
figBaseName = 'respreli_drivencells_success';
%%
cellsSubgroup1 = find(ismember(cells,intersect(cells,cellsPrefZero)));
cellsSubgroup2 = find(ismember(cells,intersect(cells,cellsPrefNinety)));
cellsSubgroup3 = find(ismember(cells,intersect(cells,setdiff(cells,cat(1,cellsPrefZero,cellsPrefNinety)))));
%% find signal correlation accross trials for each cell
trialcorrVmean = zeros(length(cycles),length(cells));
trialcorrAVmean = zeros(length(cycles),length(cells));
trialcorrVerr = zeros(length(cycles),length(cells));
trialcorrAVerr = zeros(length(cycles),length(cells));
for icyc = 1:length(cycles)
    trialcorrV{icyc} = [];
    trialcorrAV{icyc} = [];
    data = cycDataDFoverF_cmlvNoTarget{icyc};
    V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    temptrialcorrV = zeros(length(V_cycInd),length(V_cycInd),length(cells));
    temptrialcorrAV = zeros(length(AV_cycInd),length(AV_cycInd),length(cells));
    for icell = 1:length(cells)
        temptrialcorrV(:,:,icell) = corr(squeeze(data(:,cells(icell),V_cycInd)));
        temptrialcorrV(:,:,icell) = tril(temptrialcorrV(:,:,icell),-1);
        Vvector = temptrialcorrV(:,:,icell);
        Vvector = Vvector(:);
        trialcorrVmean(icyc,icell) = mean(Vvector(Vvector ~=0));
        trialcorrVerr(icyc,icell) = (std(Vvector(Vvector ~=0)))/(sqrt(length(Vvector(Vvector ~=0))));
        temptrialcorrAV(:,:,icell) = corr(squeeze(data(:,cells(icell),AV_cycInd)));
        temptrialcorrAV(:,:,icell) = tril(temptrialcorrAV(:,:,icell),-1);
        AVvector = temptrialcorrAV(:,:,icell);
        AVvector = AVvector(:);
        trialcorrAVmean(icyc,icell) = mean(AVvector(AVvector ~=0));
        trialcorrAVerr(icyc,icell) = (std(AVvector(AVvector ~=0)))/(sqrt(length(AVvector(AVvector ~=0))));
    end
    trialcorrV{icyc} = temptrialcorrV;
    trialcorrAV{icyc} = temptrialcorrAV;
end

%% scatter signal corr for aud vs vis trials
figure;
for icyc = 1:length(cycles)
    V_cycInd = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    AV_cycInd = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    subplot(3,4,icyc)
    ploterr(trialcorrVmean(icyc,cellsSubgroup1),trialcorrAVmean(icyc,cellsSubgroup1),trialcorrVerr(icyc,cellsSubgroup1),trialcorrAVerr(icyc,cellsSubgroup1),'go');
    hold on
    ploterr(trialcorrVmean(icyc,cellsSubgroup2),trialcorrAVmean(icyc,cellsSubgroup2),trialcorrVerr(icyc,cellsSubgroup2),trialcorrAVerr(icyc,cellsSubgroup2),'co');
    hold on
    ploterr(trialcorrVmean(icyc,cellsSubgroup3),trialcorrAVmean(icyc,cellsSubgroup3),trialcorrVerr(icyc,cellsSubgroup3),trialcorrAVerr(icyc,cellsSubgroup3),'ko');
    hold on
    scatter(mean(trialcorrVmean(icyc,cellsSubgroup1)),mean(trialcorrAVmean(icyc,cellsSubgroup1)),'g','filled');
    hold on
    scatter(mean(trialcorrVmean(icyc,cellsSubgroup2)),mean(trialcorrAVmean(icyc,cellsSubgroup2)),'c','filled');
    hold on
    scatter(mean(trialcorrVmean(icyc,cellsSubgroup3)),mean(trialcorrAVmean(icyc,cellsSubgroup3)),'k','filled');
	hold on
    refline(1,0);
    xlabel('visual')
    ylabel('auditory')
    title({[num2str(length(V_cycInd)) 'visual & ']; [num2str(length(AV_cycInd)) 'aud trials; ' num2str(icyc) ' cyc']});
    xlim([-0.1 0.8])
    ylim([-0.1 0.8])
    axis('square')
end

suptitle([mouse '; ' date '; ' cellgroupname])

print([fnout ['\' figBaseName '_singalcorrscatter' '.pdf']], '-dpdf')

%%

