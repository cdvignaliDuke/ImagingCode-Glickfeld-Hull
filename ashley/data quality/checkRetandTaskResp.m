%% run this after doing direction tuning analysis (need cell mask)
awFSAVdatasets_temp
% iexp = 4;
for iexp = [4 5]%1:length(expt)
% try
SubNum = expt(iexp).SubNum;
date = expt(iexp).date;
mouse = expt(iexp).mouse;

taskTimeAll = expt(iexp).time_mat;
taskFolderAll = expt(iexp).runs;

taskTime = expt(iexp).time_mat(1,:);
taskFolder = expt(iexp).runs(1,:);
taskFName = [expt(iexp).runs(1,:) '_000_000'];

retTime = expt(iexp).rettuning{2,:};
retFolder = expt(iexp).rettuning{1,:};
retFName = [expt(iexp).rettuning{1,:} '_000_000'];

% SubNum = '614';
% date = '150701';
% mouse = 'AW14';
% 
% taskTime = '1556';
% taskFolder = '002';
% taskFName = '002_000_000';
% 
% retTime = '1650';
% retFolder = '005';
% retFName = '005_000_000';

% dirTime = '1659';
% dirFolder = '005';

dataAvgFrames = 401:500;

% taskRespCutoff = 0.025;
% retRespCutoff = 0.3;
exptSummary = figure;
FOV_taskDrivenResp
FOV_retDrivenResp
quickBehaviorSummary

figure(exptSummary)
suptitle([mouse '-' date '; ' posStr '; ' sizeStr])
set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'exptQualChkSummary'), '-dpdf');
% catch
%     continue
% end
end
% 

