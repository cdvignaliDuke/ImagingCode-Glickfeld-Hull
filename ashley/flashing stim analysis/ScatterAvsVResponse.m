edit FlashingStim_dataSortedByCycle_combineDatasets.m

for icyc = 1:length(cycles)
    trialOutcome_indAll = [];
    L = ceil(30+ (cycles(icyc))*cycTime);
    cyc_ind = icyc:length(cycles);
    for i = cyc_ind
        trials = find(tCyclesOn == cycles(i));
        trialOutcome_ind = trialOutcome(trials);
        trialOutcome_indAll = cat(2,trialOutcome_indAll,trialOutcome_ind);
    end        
    cycTrialOutcome{icyc} = trialOutcome_indAll;
end

%% scatter plot

figure;
start = 1;
for icyc = 4:length(cycles)
    data = cycDataDFoverF_cmlvNoTarget{icyc};
%     V_ind = intersect(cycV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
%     AV_ind = intersect(cycAV_ind{icyc},find(strcmp(cycTrialOutcome{icyc},'success')));
    V_ind = cycV_ind{icyc};
    AV_ind = cycAV_ind{icyc};
    V_avg = squeeze(mean(mean(data(end-cycTime:end,:,V_ind),3),1));
    AV_avg = squeeze(mean(mean(data(end-cycTime:end,:,AV_ind),3),1));
    subplot(2,3,start)
    scatter(V_avg,AV_avg);
    refline(1,0);
    xlabel('visual')
    ylabel('auditory')
    title([num2str(length(V_ind)) 'visual & ' num2str(length(AV_ind)) ' auditory trials']);
    start = start+1;
end

%% color-code cells by orientation preference
DirFolder = '007';
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

%%
% set cell types, length of trial
cellGroup = baselineStimRespIndex_V;
cellsSubgroup1 = find(ismember(cellGroup,intersect(cellGroup,cellsPrefZero)));
cellsSubgroup2 = find(ismember(cellGroup,intersect(cellGroup,cellsPrefNinety)));
cellsSubgroup3 = find(ismember(cellGroup,intersect(cellGroup,setdiff([1:nCells],cat(1,cellsPrefZero,cellsPrefNinety)))));


% scatter plot
figure;
start = 1;
for icyc = 4:length(cycles)
    data = cycDataDFoverF_cmlvNoTarget{icyc};
    

    data = data(end-29:end,:,:); %%%%% CHANGE THIS TO SELECT FRAMES TO ANALYZE
    
    V_ind = cycV_ind{icyc};
    AV_ind = cycAV_ind{icyc};
    V_avg = squeeze(mean(mean(data(:,cellGroup,V_ind),3),1));
    AV_avg = squeeze(mean(mean(data(:,cellGroup,AV_ind),3),1));
    errbar_V = (std(squeeze(mean(data(:,cellGroup,V_ind),1)),[],2))./(sqrt(size(data(:,cellGroup,V_ind),3)));
    errbar_AV = (std(squeeze(mean(data(:,cellGroup,AV_ind),1)),[],2))./(sqrt(size(data(:,cellGroup,AV_ind),3)));
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
    title([num2str(length(V_ind)) 'visual & ' num2str(length(AV_ind)) ' auditory trials; ' num2str(icyc) ' cyc']);
    xlim([-0.1 0.3])
    ylim([-0.1 0.3])
    axis('square')
    start = start+1;
end

%% find modulated cells
figure;
start = 1;
for icyc = 4:length(cycles)
    data = cycDataDFoverF_cmlvNoTarget{icyc};
    

    data = data(end-29:end,:,:); %%%%% CHANGE THIS TO SELECT FRAMES TO ANALYZE
    
    V_ind = cycV_ind{icyc};
    AV_ind = cycAV_ind{icyc};
    V_avg = squeeze(mean(mean(data(:,cellGroup,V_ind),3),1));
    AV_avg = squeeze(mean(mean(data(:,cellGroup,AV_ind),3),1));
    
    
    ttestTrialsEnd = min([length(V_ind) length(AV_ind)]);
    modulationTtest = ttest(squeeze(mean(data(:,cellGroup,V_ind(1:ttestTrialsEnd)),1)),squeeze(mean(data(:,cellGroup,AV_ind(1:ttestTrialsEnd)),1)),'alpha', 0.1,'dim',2);
    modulatedCells{start} = cellGroup(find(modulationTtest == 1));
    
    errbar_V = (std(squeeze(mean(data(:,cellGroup,V_ind),1)),[],2))./(sqrt(size(data(:,cellGroup,V_ind),3)));
    errbar_AV = (std(squeeze(mean(data(:,cellGroup,AV_ind),1)),[],2))./(sqrt(size(data(:,cellGroup,AV_ind),3)));
    
    subplot(2,4,start)
    
    ploterr(V_avg(:,intersect(find(modulationTtest == 1),cellsSubgroup1)),AV_avg(:,intersect(find(modulationTtest == 1),cellsSubgroup1)),errbar_V(intersect(find(modulationTtest == 1),cellsSubgroup1),:),errbar_AV(intersect(find(modulationTtest == 1),cellsSubgroup1),:),'go');
    hold on
    ploterr(V_avg(:,intersect(find(modulationTtest == 1),cellsSubgroup2)),AV_avg(:,intersect(find(modulationTtest == 1),cellsSubgroup2)),errbar_V(intersect(find(modulationTtest == 1),cellsSubgroup2),:),errbar_AV(intersect(find(modulationTtest == 1),cellsSubgroup2),:),'co');
    hold on
    ploterr(V_avg(:,intersect(find(modulationTtest == 1),cellsSubgroup3)),AV_avg(:,intersect(find(modulationTtest == 1),cellsSubgroup3)),errbar_V(intersect(find(modulationTtest == 1),cellsSubgroup3),:),errbar_AV(intersect(find(modulationTtest == 1),cellsSubgroup3),:),'ko');
    hold on
    
    refline(1,0);
    xlabel('visual')
    ylabel('auditory')
    title([num2str(length(V_ind)) 'visual & ' num2str(length(AV_ind)) ' auditory trials']);
    xlim([-0.1 0.3])
    ylim([-0.1 0.3])
    axis('square')
    start = start+1;
end

