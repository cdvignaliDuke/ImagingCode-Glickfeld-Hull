function [cellsSelect, OSI, DSI] = OriCellSets(rc, expt, iexp,cellsOnly);
%display experiment
disp([num2str(expt(iexp).date) ' i' num2str(expt(iexp).SubNum)])
%load direction tuning data
dataPath = fullfile(rc.ashleyAnalysis,expt(iexp).mouse,expt(iexp).folder, expt(iexp).date, expt(iexp).dirtuning);
if cellsOnly == 1
    load(fullfile(dataPath,'cell&dendriteTC.mat'))
    data_tc_subnp = cell_tc;
elseif cellsOnly == 2    
    load(fullfile(dataPath,'timecourses_dendrites.mat'))
else
    load(fullfile(dataPath, 'timecourses.mat'));
end

%load mworks file
mworks = ['data-' 'i' expt(iexp).SubNum '-' expt(iexp).date '-' expt(iexp).dirtuning_time]; 
load (fullfile(rc.pathStr,mworks));

%params
nON = (input.nScansOn)./10;
nOFF = (input.nScansOff)./10;
nStim = input.gratingDirectionStepN;
nTrials = input.trialSinceReset;
% nRep = size(data_tc_subnp,1)./((nON+nOFF)*nStim);
% nTrials = (nStim.*nRep);
data_tc_subnp = data_tc_subnp(1:(nON+nOFF)*nTrials,:);
DirectionDeg = cell2mat(input.tGratingDirectionDeg);
Dirs = unique(DirectionDeg);
base_win = 6:10;
resp_win = 12:16;



%find dF/F of responses for all cells
stimOFF_ind = 1:nOFF+nON:size(data_tc_subnp,1);
dF_data = zeros(size(data_tc_subnp));
dFoverF_data = zeros(size(data_tc_subnp));
for i = 1:nTrials
    indAll = stimOFF_ind(i):stimOFF_ind(i)+(nOFF+nON-1);
    indF = stimOFF_ind(i)+5:stimOFF_ind(i)+(nOFF-1);
    dF_data(indAll,:) = bsxfun(@minus,data_tc_subnp(indAll,:),mean(data_tc_subnp(indF,:),1));
    dFoverF_data(indAll,:) = bsxfun(@rdivide,dF_data(indAll,:),mean(data_tc_subnp(indF,:),1));
end

%divide by trials
stimON_ind = nOFF+1:nOFF+nON:size(dFoverF_data,1);
dFoverFCellsTrials = zeros(10+nON,size(dFoverF_data,2),nTrials);
for i = 1:nTrials
    dFoverFCellsTrials(:,:,i) = dFoverF_data(stimON_ind(i)-10:stimON_ind(i)+(nON-1),:);
end

%average direction tuning curves and find sig resp cells
dFoverF_DirResp_avg = zeros(nStim, size(dFoverFCellsTrials,2));
dFoverF_DirResp_sem = zeros(nStim, size(dFoverFCellsTrials,2));
h_dir = zeros(nStim, size(dFoverFCellsTrials,2));
if size(dFoverFCellsTrials,3) > length(DirectionDeg)
   nTrials = length(DirectionDeg); 
end

for i = 1:nStim
    trials = find(DirectionDeg(:,1:nTrials) == Dirs(i));
    dFoverF_DirResp_avg(i,:) = squeeze(mean(mean(dFoverFCellsTrials(resp_win,:,trials),1),3) - mean(mean(dFoverFCellsTrials(base_win,:,trials),1),3));
    dFoverF_DirResp_sem(i,:) = squeeze(std(mean(dFoverFCellsTrials(resp_win,:,trials),1) - mean(dFoverFCellsTrials(base_win,:,trials),1), [],3))./(sqrt(length(trials)));
    [h_dir(i,:), p] = ttest(mean(dFoverFCellsTrials(resp_win,:,trials),1), mean(dFoverFCellsTrials(base_win,:,trials),1), 'dim', 3, 'tail', 'right', 'alpha', 0.05); 
end
resp_dir_ind = find(sum(h_dir,1)>0);

%average orientation tuning curves and find sig resp cells
h_ori = zeros(nStim/2, size(dFoverFCellsTrials,2));
dFoverF_OriResp_avg = zeros(nStim/2, size(dFoverFCellsTrials,2));
dFoverF_OriResp_sem = zeros(nStim/2, size(dFoverFCellsTrials,2));
dFoverF_OriResp_TC = zeros(nStim/2,size(dFoverFCellsTrials,1),size(dFoverFCellsTrials,2));
for i = 1:nStim/2
    trials = find(DirectionDeg(:,1:nTrials) == Dirs(i) | DirectionDeg(:,1:nTrials) == Dirs(i+nStim/2));
    dFoverF_OriResp_avg(i,:) = squeeze(mean(mean(dFoverFCellsTrials(resp_win,:,trials),1),3) - mean(mean(dFoverFCellsTrials(base_win,:,trials),1),3));
    dFoverF_OriResp_TC(i,:,:) = mean(dFoverFCellsTrials(:,:,trials),3);
    dFoverF_OriResp_sem(i,:) = squeeze(std(mean(dFoverFCellsTrials(resp_win,:,trials),1) - mean(dFoverFCellsTrials(base_win,:,trials),1), [],3))./(sqrt(length(trials)));
    [h_ori(i,:), p] = ttest(mean(dFoverFCellsTrials(resp_win,:,trials),1), mean(dFoverFCellsTrials(base_win,:,trials),1), 'dim', 3, 'tail', 'right', 'alpha', 0.05); 
end
resp_ori_ind = find(sum(h_ori,1)>0);
resp_ind = unique([resp_dir_ind resp_ori_ind]);

%bootstrap to find reliable tuning
dFoverF_DirResp_avg_boot = zeros(nStim, size(dFoverFCellsTrials,2), 1000);
dFoverF_OriResp_avg_boot = zeros(nStim/2, size(dFoverFCellsTrials,2), 1000);
for i = 1:nStim
    trials = find(DirectionDeg(:,1:nTrials) == Dirs(i));
    if i <= nStim/2
        trials_ori = find(DirectionDeg(:,1:nTrials) == Dirs(i) | DirectionDeg(:,1:nTrials) == Dirs(i+nStim/2));
    end
    for ii = 1:1000;
        samp = randsample(trials, length(trials), 1);
        dFoverF_DirResp_avg_boot(i,:,ii) = squeeze(mean(mean(dFoverFCellsTrials(resp_win,:,samp),1),3) - mean(mean(dFoverFCellsTrials(base_win,:,samp),1),3));
        if i <= nStim/2
            samp = randsample(trials_ori, length(trials_ori), 1);
            dFoverF_OriResp_avg_boot(i,:,ii) = squeeze(mean(mean(dFoverFCellsTrials(resp_win,:,samp),1),3) - mean(mean(dFoverFCellsTrials(base_win,:,samp),1),3));
        end
    end
end
[boot_dir_max boot_dir_ind] = max(dFoverF_DirResp_avg_boot,[],1);
[boot_ori_max boot_ori_ind] = max(dFoverF_OriResp_avg_boot,[],1);

%find those cells with >90% CI of having Ori preference 
n_dir = zeros(nStim, size(dFoverFCellsTrials,2));
n_ori = zeros(nStim, size(dFoverFCellsTrials,2));
ind_rank = zeros(1, size(dFoverFCellsTrials,2));
max_dir_ind = nan(1, size(dFoverFCellsTrials,2));
max_ori_ind = nan(1, size(dFoverFCellsTrials,2));
no_ori_ind = nan(1, size(dFoverFCellsTrials,2));
for i= 1:size(dFoverFCellsTrials,2)
    [n_dir(:,i)] = histc(squeeze(boot_dir_ind(1,i,:)),[1:nStim]);
    [n_ori(:,i)] = histc(squeeze(boot_ori_ind(1,i,:)),[1:nStim]);
    [a, b] =  sort(n_dir(:,i), 'descend');
    if a(1) > 900 
    	max_dir_ind(:,i) = b(1);
    elseif abs(b(1)-b(2)) == nStim/2  & a(1)+a(2)>900
        max_dir_ind(:,i) = b(1);
    end
    [a, b] =  sort(n_ori(:,i), 'descend');
    if a(1) > 900
        max_ori_ind(:,i) = b(1);
    end
    if a(nStim/2) >= 100
        no_ori_ind(:,i) = 1;
    end
end

%find cells with <10% CI of having Ori preference, but still driven
untuned_ind = intersect(find(~isnan(no_ori_ind)),resp_ind);

%determine OSI/DSI for reliable cells
x(1,:) = [1:nStim];
x(2,:) = [3 4 1 2 7 8 5 6];
x(3,:) = [5 6 7 8 1 2 3 4];
pref_val = nan(1, size(dFoverFCellsTrials,2));
ori_null_val = nan(1, size(dFoverFCellsTrials,2));
dir_null_val = nan(1, size(dFoverFCellsTrials,2));
dFoverF_DirResp_avg_rect = dFoverF_DirResp_avg;
dFoverF_DirResp_sem_rect = dFoverF_DirResp_sem;
dFoverF_DirResp_avg_rect(find(dFoverF_DirResp_avg<0))= 0;
dFoverF_DirResp_sem_rect(find(dFoverF_DirResp_avg<0))= NaN;
dFoverF_OriResp_avg_rect = dFoverF_OriResp_avg;
dFoverF_OriResp_sem_rect = dFoverF_OriResp_sem;
dFoverF_OriResp_avg_rect(find(dFoverF_OriResp_avg<0))= 0;
dFoverF_OriResp_sem_rect(find(dFoverF_OriResp_avg<0))= NaN;
for i = 1:size(dFoverFCellsTrials,2)
    if ~isnan(max_dir_ind(:,i))
        pref_val(i) = dFoverF_DirResp_avg_rect(max_dir_ind(1,i),i);
        ori_null_val(i) = dFoverF_DirResp_avg_rect(x(2,max_dir_ind(1,i)),i);
        dir_null_val(i) = dFoverF_DirResp_avg_rect(x(3,max_dir_ind(1,i)),i);
    elseif ~isnan(max_ori_ind(:,i))
        pref_val(i) = dFoverF_OriResp_avg_rect(max_ori_ind(1,i),i);
        ori_null_val(i) = dFoverF_OriResp_avg_rect(x(2,max_ori_ind(1,i)),i);
    end
end
OSI = (pref_val-ori_null_val)./(pref_val+ori_null_val);
DSI = (pref_val-dir_null_val)./(pref_val+dir_null_val);

%find preferred ori for dir selective cells
dir_ori_ind = max_dir_ind;
dir_ori_ind(find(max_dir_ind>4)) = max_dir_ind(find(max_dir_ind>4))-4;
ori_ind_all = dir_ori_ind;
ori_ind_all(find(isnan(dir_ori_ind))) = max_ori_ind(find(isnan(dir_ori_ind)));

%create cell array for ori selective cells
for iOri = 1:nStim/2
    cellsSelect{iOri} = intersect(find(ori_ind_all == iOri), find(OSI>0.3));
end

cellsSelect{iOri+1} = untuned_ind;
cellsSelect{iOri+2} = intersect(find(~isnan(ori_ind_all)), find(OSI>0.3));

figure;
h = histogram(OSI,10);
xlabel('OSI')
ylabel('n cells')
title([expt(iexp).mouse '-' expt(iexp).date])
print(fullfile(dataPath,'OSI_hist'),'-dpdf')

if cellsOnly == 1
    save(fullfile(dataPath, 'cellsSelect_cellsOnly.mat'), 'cellsSelect', 'OSI', 'DSI','ori_ind_all','max_dir_ind','dFoverF_OriResp_avg_rect','dFoverF_OriResp_sem_rect','dFoverF_DirResp_avg_rect','dFoverF_DirResp_sem_rect','dFoverF_OriResp_TC');
elseif cellsOnly == 2
    save(fullfile(dataPath, 'cellsSelect_dendrites.mat'), 'cellsSelect', 'OSI', 'DSI','ori_ind_all','max_dir_ind','dFoverF_OriResp_avg_rect','dFoverF_OriResp_sem_rect','dFoverF_DirResp_avg_rect','dFoverF_DirResp_sem_rect','dFoverF_OriResp_TC');
else
    save(fullfile(dataPath, 'cellsSelect.mat'), 'cellsSelect', 'OSI', 'DSI','ori_ind_all','max_dir_ind','dFoverF_OriResp_avg_rect','dFoverF_OriResp_sem_rect','dFoverF_DirResp_avg_rect','dFoverF_DirResp_sem_rect','dFoverF_OriResp_TC');
end

end






