%% run this after doing direction tuning analysis (need cell mask)
awFSAVdatasets
for iexp = 1:length(expt)
% try
SubNum = expt(iexp).SubNum;
date = expt(iexp).date;
mouse = expt(iexp).mouse;

taskTimeAll = expt(iexp).time_mat;
taskFolderAll = expt(iexp).runs;

taskTime = expt(iexp).time_mat(2,:);
taskFolder = expt(iexp).runs(2,:);
taskFName = [expt(iexp).runs(2,:) '_000_000'];

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

FOV_taskDrivenResp
FOV_retDrivenResp
quickBehaviorSummary
toc
% catch
%     continue
% end
end


