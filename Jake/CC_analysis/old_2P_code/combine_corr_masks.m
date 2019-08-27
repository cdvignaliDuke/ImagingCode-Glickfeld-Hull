function sm_out = combine_corr_masks(img_reg, pixel_mask, threshold)
%function for combining IC masks based upon their activity correlations.
%Written for the CRP analysis

%allocate memory and define dimensions
[npw, nph, ~] = size(img_reg);
this_mask = 1;
if length(size(pixel_mask))==2 %break up the 2D mask into a 3D mask
    nmask = unique(pixel_mask);
    nmask = nmask(nmask >0);
    sm = zeros(npw,nph,length(nmask));
    for ii=nmask %1:nmask
        sm_placeholder = zeros(npw,nph);
        sm_placeholder(pixel_mask==ii) = 1;
        sm(:,:,this_mask)= sm_placeholder;
        this_mask = this_mask +1;
    end
elseif length(size(pixel_mask))==3
    nmask = [1:size(pixel_mask,3)];
    sm = pixel_mask;
end

%get TCs from img_reg
tc_avg = getTC(img_reg, sm, length(nmask));

%run correlation on the TCs. 
data_corr = corrcoef(tc_avg);

%consolidate timecourses that are highly correlated
[i, j] = find(and(data_corr>threshold,data_corr<1)); %look for neurons which have sig corr but are not on the unity line
sm_temp = reshape(sm,npw*nph,length(nmask));
remove_2 =[];
if ~isempty(i)
    comp = [j(1); i(j == j(1))];
    total = unique(i); %its possible for a neuron to have sig corr with more than one other neuron. So find the total unique neurons with sig corrs.
    rest_ind = ismember(total, comp); %find all the masks not included in the comparison group
    rest_comp = total(~rest_ind); %
    
    %add all the masks in the comparison group to sm(:,:,comp(1))
    for k = 2:length(comp)
        sm_temp(find(sm(:,:,comp(k))),comp(1)) = 1; %put all of the masks from the comp group onto one mask
        sm(:, :, comp(1)) = reshape(sm_temp(:,comp(1)), npw ,nph); %update sm
    end
    
    if length(rest_comp) == 2
        sm_temp( find(sm(:,:,rest_comp(2))) ,rest_comp(1) ) = 1; %also add one of the rest comp masks to the all comp mask 
        sm(:, :, rest_comp(1)) = reshape(sm_temp(:,rest_comp(1)), npw ,nph); %update that sm mask
        remove_2 = rest_comp(2);
    elseif length(rest_comp) > 2
        for k = 1:length(rest_comp)
            if k > length(rest_comp)
                break
            else
                rest_comp1 = rest_comp(k);
                rest_comp2 = j(i == rest_comp1); %find all the masks which sig corr with rest_comp(k)
                rest_comp(ismember(rest_comp,rest_comp2))=[]; %removes rest_comp2 from rest_comp
                for kk = 1:length(rest_comp2)
                    sm_temp(find(sm(:,:,rest_comp2(kk))),rest_comp1) = 1; %put all of the masks from the comp group onto one mask
                    sm(:, :, rest_comp1) = reshape(sm_temp(:,rest_comp1), npw ,nph);
                end
                remove_2 = [remove_2; rest_comp2]; %keep track of 
            end
        end
    end
    remove_2 = [comp(2:end);remove_2]; 
    sm(:, :, remove_2) = [];
end

%% combine masks which are highly overlapping

%sort mask_cell according to mask size
sm = reshape(sm, npw*nph, size(sm,3));
sm_sums = sum(sm,1);
[~,sorted_ind] = sort(sm_sums, 'descend');
sm = sm(:,[sorted_ind]);

%finds overlapping pixels of ROIs and based on amount of overlap decides whether
%to group them or to split them- if splitting, then overlapping pixels are
%eliminated from both ROIs
mask_cell_temp = sm;
nIC = size(mask_cell_temp,2);

for ic = 1:nIC
    %make a mask which includes all masks except the current one
    comp_mask = mask_cell_temp;
    comp_mask(:,ic) = 0;
    
    %find which masks overlap with current mask
    overlap_mask = comp_mask + repmat(mask_cell_temp(:,ic),1,nIC);
    overlap_ICs = max(overlap_mask,[],1);
    overlap_ICs = find(overlap_ICs>1);
    if isempty(overlap_ICs)
        continue
    end
    
    %combine the highly overlapping masks
    for mask_num = overlap_ICs
        this_mask_size = sum(mask_cell_temp(:,ic),1);
        overlap_ratio = length(find(overlap_mask(:,mask_num)==2))/this_mask_size;
        if overlap_ratio > 0.9
            mask_cell_temp(:,ic) = logical(mask_cell_temp(:,ic) + mask_cell_temp(:,mask_num));
            mask_cell_temp(:,mask_num) = 0;
        end
    end
end
sm_out = reshape(mask_cell_temp, npw, nph, size(mask_cell_temp,2));





