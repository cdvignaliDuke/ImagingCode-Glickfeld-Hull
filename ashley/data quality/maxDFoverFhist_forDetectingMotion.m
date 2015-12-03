ialign = 1;
awFSAVdatasets
for iexp = 1:length(expt)
% iexp = 9;
    divideupdatabyalignment
cycDataDFoverF_meanCells = cellfun(@squeeze,cellfun(@mean,cycDataDFoverF,num2cell(repmat(2,1,length(cycDataDFoverF))),'UniformOutput',0),'UniformOutput',0);
cycDataDFoverF_absDiffPerTrial = cellfun(@abs,cellfun(@diff,cycDataDFoverF_meanCells,num2cell(repmat(1,1,length(cycDataDFoverF_meanCells))),num2cell(repmat(1,1,length(cycDataDFoverF_meanCells))),'UniformOutput',0),'UniformOutput',0);
cycDataDFoverF_maxDiffPerTrial = cellfun(@max,cycDataDFoverF_absDiffPerTrial,cell(1,length(cycDataDFoverF_absDiffPerTrial)),num2cell(repmat(1,1,length(cycDataDFoverF_meanCells))),'UniformOutput',0);

% for i = 1:length(cycDataDFoverF_maxDiffPerTrial)
%     subplot(3,4,i)
%     hist(cycDataDFoverF_maxDiffPerTrial{i},100)
% end

% bmDist = gmdistribution.fit([cycDataDFoverF_maxDiffPerTrial{:}]',2);
% highDistMean = bmDist.mu(2);
% highDistSTD = squeeze(bmDist.Sigma(:,:,2));
% motionThreshold = highDistMean-(2*highDistSTD);

figure; hist([cycDataDFoverF_maxDiffPerTrial{:}],100)
hold on
% vline(motionThreshold,'k')
title([expt(iexp).mouse '-' expt(iexp).date '; max dF/F'])

set(0,'defaultfigurepaperorientation','portrait');
set(0,'defaultfigurepapersize',[8.5 11]);
set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);
fnout = ['Z:\Analysis\' mouse '\' expFolder '\' dateofexp '\' runstr];
try 
    cd([fnout '\maxDFoverF_hist'])
catch
    mkdir(fnout,'maxDFoverF_hist')
    cd([fnout '\maxDFoverF_hist'])
end
% print([fnout '\maxDFoverF_alltrials_hist'],'-dpdf')

%%
thrsh = 0.02;
cycDataDFoverF_indTrialsMaxOverThreshold = cellfun(@(y)find(y),cellfun(@(x) x > thrsh,cycDataDFoverF_maxDiffPerTrial,'UniformOutput',0),'UniformOutput',0);
save('probMotionTrials','cycDataDFoverF_indTrialsMaxOverThreshold');
%%
end