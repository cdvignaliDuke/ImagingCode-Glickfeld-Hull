fileSave = fullfile(rc.ashleyAnalysis,mouse_name,folder, date_name, dir_run);
cd(fileSave);
load('TuningPreferences.mat')

% dataTrialStart = cycDataDFoverF_cmlvNoTarget{4};
% v_ind = cycV_ind{4};
% % a_ind = cycA_ind{1};
% a_ind = cycAV_ind{4};
% 
% drivencells = find(any(mean(cycDataDFoverF_cmlvNoTarget{end},3) >0.05,1));
% 
% preStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
% for itrial =1:size(v_ind,1);
%     for icell = 1:size(dataTrialStart,2)
%         preStimResp_V(itrial,icell) = mean(dataTrialStart(1:30,icell,v_ind(itrial)),1);
%     end
% end
% 
% baselineStimResp_V = zeros(size(v_ind,2),size(dataTrialStart,2));
% for itrial = 1:size(v_ind,1);
%     for icell = 1:size(dataTrialStart,2)
%         baselineStimResp_V(itrial,icell) = mean(dataTrialStart(36:end,icell,v_ind(itrial)),1);
%     end
% end
% 
% baselineStimRespTtest_V= ttest(preStimResp_V,baselineStimResp_V,'alpha', 0.01);
% baselineStimRespIndex_V = find(baselineStimRespTtest_V == 1);
nOri = length(unique(oriPref_ind));
Oris = 0:180/nOri:180-(180/nOri);
for iOri = 1:nOri
cellsPref{iOri} = find(dirPref_ind == iOri | dirPref_ind == iOri+4);
cellsSelect_dir{iOri} = intersect(dirSlctvCells,cellsPref{iOri});
cellsSelect_ori{iOri} = intersect(oriSlctvCells,cellsPref{iOri});
cellsSelect{iOri} = union(cellsSelect_dir{iOri},cellsSelect_ori{iOri});
% cellsPrefResp{iOri} = intersect(baselineStimRespIndex_V,cellsPref{iOri});
end

nCells = size(dataTC,2);
oriSlctvCellsAll = union(oriSlctvCells,dirSlctvCells);
notSlctvCells = setdiff([1:nCells],oriSlctvCellsAll);
% notRespCells = setdiff([1:nCells],baselineStimRespIndex_V);
