%% run this after doing direction tuning analysis (need cell mask)
awFSAVdatasets_naiveRandBaseOri
iexp = 1;
% for iexp = 4:5
% % try
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

% SubNum = '625';
% date = '160301';
% mouse = 'AW25';
% 
% taskTime = '1205';
% taskFolder = '004';
% taskFName = '004_000_000';
% 
% retTime = '1205';
% retFolder = '004';
% retFName = '004_000_000';
% 
% dirTime = '1659';
% dirFolder = '005';

dataAvgFrames = 401:500;

% taskRespCutoff = 0.025;
% retRespCutoff = 0.3;

set(0,'defaultfigurepaperorientation','landscape');
set(0,'defaultfigurepapersize',[11 8.5]);
set(0,'defaultfigurepaperposition',[.25 .25 [11 8.5]-0.5]);
exptSummary = figure;
FOV_taskDrivenResp
figure(exptSummary)
suptitle([mouse '-' date '; ' posStr '; ' sizeStr])
FOV_retDrivenResp
quickBehaviorSummary

figure(exptSummary)
print(fullfile('Z:\analysis\',mouse,'two-photon imaging', date,'exptQualChkSummary'), '-dpdf');
% catch
%     continue
% end
% close all
% end
% 

