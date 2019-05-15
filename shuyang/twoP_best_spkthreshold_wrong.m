% This function decides the most optimal threshold when identifying spikes
% i.e: how many std of the 1st derivatives of raw fluorescence
% used in twoP_MDspikes, 
% output:  different ave baseline (stationary) firing rates w/ different thresholds, 
% average FR for all cells during stationary due to the best threshold
% (when ave FR is around 1).

function [avesWvarystd, best_thres,aveFR_allcells_stay,std_best] = best_spkthreshold_2P (deriv, frm_stay,std2,std2_5,std3,std_deriv) 
nspk2 = zeros(1,size(deriv,2)); % num of cells
frate2 = zeros(size(frm_stay,2),size(deriv,2)); % num of stationary windows*num of cells
for i = 1: size(frm_stay,2) %for each stationary window
    t = length(frm_stay{i})/30; %t = nframes*1/30s, thus the unit of t is second
    der_stay_temp = deriv(frm_stay{i},:); 
    for j = 1: size(der_stay_temp,2) % for each cell
        nspk2(j) = sum(der_stay_temp(:,j) >= std2(j)); %number of spikes during this window for each cell
    end
    frate2(i,:) = nspk2/t; % 2 dimensions: window index*cell
end
aveFR_cells2 = mean(frate2);%ave firing rate for all cells across all running windows (1*cells)
aveFR_allcells2 = mean(aveFR_cells2); %ave firing rate across all cells (1 single num)

nspk1 = zeros(1,size(deriv,2)); 
frate1 = zeros(size(frm_stay,2),size(deriv,2)); 
for i = 1: size(frm_stay,2) 
    t = length(frm_stay{i})/30; 
    der_stay_temp = deriv(frm_stay{i},:); 
    for j = 1: size(der_stay_temp,2) 
        nspk1(j) = sum(der_stay_temp(:,j) >= std_deriv(j)); 
    end
    frate1(i,:) = nspk1/t; 
end
aveFR_cells1 = mean(frate1);
aveFR_allcells1 = mean(aveFR_cells1); 

nspk3 = zeros(1,size(deriv,2)); 
frate3 = zeros(size(frm_stay,2),size(deriv,2)); 
for i = 1: size(frm_stay,2) 
    t = length(frm_stay{i})/30; 
    der_stay_temp = deriv(frm_stay{i},:); 
    for j = 1: size(der_stay_temp,2) 
        nspk3(j) = sum(der_stay_temp(:,j) >= std3(j)); 
    end
    frate3(i,:) = nspk3/t; 
end
aveFR_cells3 = mean(frate3);
aveFR_allcells3 = mean(aveFR_cells3); 

nspk2_5 = zeros(1,size(deriv,2)); 
frate2_5 = zeros(size(frm_stay,2),size(deriv,2)); 
for i = 1: size(frm_stay,2) 
    t = length(frm_stay{i})/30; 
    der_stay_temp = deriv(frm_stay{i},:); 
    for j = 1: size(der_stay_temp,2) 
        nspk2_5(j) = sum(der_stay_temp(:,j) >= std2_5(j)); 
    end
    frate2_5(i,:) = nspk2_5/t; 
end
aveFR_cells2_5 = mean(frate2_5);
aveFR_allcells2_5 = mean(aveFR_cells2_5); 

% find the aveFR_allcells closest to 1, 
all_aves = [aveFR_allcells1,aveFR_allcells2,aveFR_allcells3,aveFR_allcells2_5];
aveFR_allcells = [aveFR_cells1;aveFR_cells2;aveFR_cells3;aveFR_cells2_5];
stds = [std_deriv;std2;std3;std2_5];
thres = [1,2,3,2.5]; % !!!!!!!! the order of the elements in these 4 lines are super important, make sure they match!
delta = abs(all_aves - 1);
bestave = all_aves(delta == min(delta));
ind_best = find(all_aves == bestave);
best_thres = thres(ind_best); % when firing rate is around 1 during stationary (baseline), that relative threshold is the most optimal.
aveFR_allcells_stay = aveFR_allcells((ind_best),:); 
std_best = stds(best_thres,:);
avesWvarystd = struct('aveFR_allcells1',aveFR_allcells1,'aveFR_allcells2',...
    aveFR_allcells2,'aveFR_allcells3',aveFR_allcells3,'aveFR_allcells2_5',aveFR_allcells2_5);


% plot aveFR cells ----------------------------------------------------------------------------------------
figure; hist(aveFR_cells2_5); title('aveFR during stationary w std 2.5');
figure; hist(aveFR_cells2); title('aveFR during stationary w std 2')

end







