% This function decides the most optimal threshold when identifying spikes
% i.e: how many std of the 1st derivatives of raw fluorescence
% used in spk_derivative_threshold, 
% output:  different ave baseline (stationary) firing rates w/ different thresholds, 
% average FR for all cells during stationary due to the best threshold
% (when ave FR is around 1).
% spk_inx: indexes of the frames that has a Ca event
% spk_bi_cellmat: a matrix with 0s and 1s, 1 is Ca event.
% std best: standard deviation of the first derivatives with the best threshold

function [aveFRsWstds, best_thres,bestFR,std_best,spk_inx,FRstay_cells,...
    spk_bi_cellmat] = twoP_best_spkthreshold (deriv, frm_stay,TCave) 
std_deriv = std(deriv);
std2 = 2*std_deriv;
std3 = 3*std_deriv;
std2_5 = 2.5*std_deriv;
% when you count spikes, there shoouldn't be any contineous indices
% get the best threshold std for deciding spikes, ave FRs during stationary for each cell
% get the peak values of derivatives and their indices for all of the cells
% you need to get the peaks first because if you just threshold the derivatives,
% there can be a few contineous derivatives bigger then the threshold, but that's 1 peak.
peaks = {};
locs = {};
for i = 1: size(deriv,2) %for each cell
    [peaks{i},locs{i}] = findpeaks(deriv(:,i));%peaks and indices of the peaks for each cell
    peaks{i} = (peaks{i})';
    locs{i} = (locs{i})';
    %peaks and locs both have number of n elements, n is num of cells
end

%For each cell: find the indices of the peaks which is bigger than the threshold, and
%count how many of these fall in the stationary peroids to calculate spike
%rates for stationary
pk_abv1_inx = {}; pk_abv2_inx = {}; pk_abv2_5_inx = {}; pk_abv3_inx = {};
spk1 = zeros(1,size(deriv,2));spk2 = zeros(1,size(deriv,2));spk2_5 = zeros(1,size(deriv,2));spk3 = zeros(1,size(deriv,2));

for i = 1: size(deriv,2) %for each cell
    pk_abv1_inx{i} = locs{i}(peaks{i} >= std_deriv(i))+1;% these are all of the frames where there's a spike, if the threshold is 1std, frame index = derivative index +1
    %pk_abv_ind has number of n elements, n is num of cells. these indices are
    %the frames that has a spike
    pk_abv2_inx{i} = locs{i}(peaks{i} >= std2(i))+1;
    pk_abv2_5_inx{i} = locs{i}(peaks{i} >= std2_5(i))+1;
    pk_abv3_inx{i} = locs{i}(peaks{i} >= std3(i))+1;
    
    %count how many spikes are in the stationary periods: get the intersection
    %of indices of peaks above the threshold and indices of frames for stationary
    spk1(i) = length(intersect(pk_abv1_inx{i},frm_stay));
    spk2(i) = length(intersect(pk_abv2_inx{i},frm_stay));
    spk2_5(i) = length(intersect(pk_abv2_5_inx{i},frm_stay));
    spk3(i) = length(intersect(pk_abv3_inx{i},frm_stay));
end

%calculate the ave spike rates of all cells based on different thresholds during stationary
t_stay = length(frm_stay)/30;
aveFR1 = mean(spk1)/t_stay;
aveFR2 = mean(spk2)/t_stay;
aveFR2_5 = mean(spk2_5)/t_stay;
aveFR3 = mean(spk3)/t_stay;

% the FR closest to 1 is the "best", save the best threshold, firing rate
% with the best thershold
thres = [1,2,2.5,3]; % !!!!!!!! the order of the elements in these 4 lines are super important, make sure they match!
all_aveFRs = [aveFR1,aveFR2,aveFR2_5,aveFR3];
all_FR_cells = [spk1/t_stay;spk2/t_stay;spk2_5/t_stay;spk3/t_stay];
stds = [std_deriv;std2;std3;std2_5];
all_pk_abv_inx = {pk_abv1_inx,pk_abv2_inx,pk_abv2_5_inx,pk_abv3_inx};
delta = abs(all_aveFRs - 1);
bestFR = all_aveFRs(delta == min(delta));
ind_best = find(all_aveFRs == bestFR);
best_thres = thres(ind_best);
std_best = stds(best_thres,:);
spk_inx = all_pk_abv_inx{ind_best};
aveFRsWstds = struct('aveFR1',aveFR1,'aveFR2',aveFR2,'aveFR2_5',aveFR2_5,'aveFR3',aveFR3);
FRstay_cells = all_FR_cells(ind_best,:);

%write a binary matrix of frames*cells. this matrix will only contain 0s
%and 1s, 1 for spike
spk_bi_cellmat = zeros(size(TCave,1),size(TCave,2));
for i = 1:size(spk_inx,2)             %for each cell
    spk_inx_cell = spk_inx{i};
    spk_bi_cellmat(spk_inx_cell,i) = 1;
end

% plot aveFR cells ----------------------------------------------------------------------------------------
figure; subplot(1,2,1); hist(spk2_5/t_stay); title('aveFR during stationary w std 2.5');
subplot(1,2,2);hist(spk2/t_stay); title('aveFR during stationary w std 2');
end


