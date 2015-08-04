DirFolder = '005';
fileSave = fullfile('Z:\analysis\',mouse,'two-photon imaging', date, DirFolder);
cd(fileSave);
% load('cellSelectivityIndices.mat')
load('TuningPreferences.mat')

cellsPrefZero = find(dirPref_ind == 1 | dirPref_ind == 5);
cellsPrefNinety = find(dirPref_ind == 3 | dirPref_ind == 7);
otherCells = setdiff([1:size(data,2)]',cat(1,cellsPrefZero,cellsPrefNinety));

figure;
start = 1;
for icyc = 4:length(cycles)
    data = cycDataDFoverF_cmlvNoTarget{icyc};
    V_ind = cycV_ind{icyc};
    AV_ind = cycAV_ind{icyc};
    V_std = squeeze(mean(std(data(:,:,V_ind),0,3),1));
    AV_std = squeeze(mean(std(data(:,:,AV_ind),0,3),1));   
    subplot(2,3,start)
    scatter(V_std(:,cellsPrefZero),AV_std(:,cellsPrefZero),'g');
    hold on
%     scatter(mean(V_std(:,cellsPrefZero)),mean(V_std(:,cellsPrefZero)),'g','filled');
%     hold on
    scatter(V_std(:,cellsPrefNinety),AV_std(:,cellsPrefNinety),'c');
    hold on
%     scatter(mean(V_std(:,cellsPrefNinety)),mean(V_std(:,cellsPrefNinety)),'c','filled');
%     hold on
    scatter(V_std(:,otherCells),AV_std(:,otherCells),'k');
    hold on
%     scatter(mean(V_std(:,otherCells)),mean(V_std(:,otherCells)),'k','filled');
%     hold on
    refline(1,0);
    xlabel('visual')
    ylabel('auditory')
    title([num2str(length(V_ind)) 'visual & ' num2str(length(AV_ind)) ' auditory trials']);
%     xlim([-0.1 0.3])
%     ylim([-0.1 0.3])
    xlim([0 0.15])
    ylim([0 0.15])
    axis('square')
    start = start+1;
end

%% responsive cells
figure;
start = 1;
for icyc = 4:length(cycles)
    data = cycDataDFoverF_cmlvNoTarget{icyc};
    V_ind = cycV_ind{icyc};
    AV_ind = cycAV_ind{icyc};
    V_std = squeeze(mean(std(data(:,:,V_ind),0,3),1));
    AV_std = squeeze(mean(std(data(:,:,AV_ind),0,3),1));   
    subplot(2,3,start)
    scatter(V_std(:,intersect(baselineStimRespIndex_V,cellsPrefZero)),AV_std(:,intersect(baselineStimRespIndex_V,cellsPrefZero)),'g');
    hold on
    scatter(mean(V_std(:,cellsPrefZero)),mean(AV_std(:,cellsPrefZero)),'g','filled');
    hold on
    scatter(V_std(:,intersect(baselineStimRespIndex_V,cellsPrefNinety)),AV_std(:,intersect(baselineStimRespIndex_V,cellsPrefNinety)),'c');
    hold on
    scatter(mean(V_std(:,cellsPrefNinety)),mean(AV_std(:,cellsPrefNinety)),'c','filled');
    hold on
    scatter(V_std(:,intersect(baselineStimRespIndex_V,otherCells)),AV_std(:,intersect(baselineStimRespIndex_V,otherCells)),'k');
    hold on
    scatter(mean(V_std(:,otherCells)),mean(AV_std(:,otherCells)),'k','filled');
    hold on
    refline(1,0);
    xlabel('visual')
    ylabel('auditory')
    title([num2str(length(V_ind)) 'visual & ' num2str(length(AV_ind)) ' auditory trials']);
%     xlim([-0.1 0.3])
%     ylim([-0.1 0.3])
    xlim([0 0.15])
    ylim([0 0.15])
    axis('square')
    start = start+1;
end