function mask_final = processMask_CRP(mask_cell)
%get dimensions of mask, allocate memory
[npw, nph, nIC] = size(mask_cell);
mask_cell_temp = zeros(npw*nph, nIC);

%identifies masks which have more than one object that dont touch.
%Seperates those objects into seperate masks. 
for ic = 1:nIC
    if length(unique(reshape(mask_cell(:,:,ic),[1 npw*nph])))>2 %if there is more than one object in this mask..
        %         data_tc_temp = stackGetTimeCourses(img_reg,mask_cell(:,:,ic));
        %         data_corr_temp = corrcoef(data_tc_temp);
        ind_rem = 1:length(unique(reshape(mask_cell(:,:,ic),[1 npw*nph])))-1;
        for i = ind_rem %for each unique object in this mask
            ind_new = find(mask_cell(:,:,ic) == i);
            if length(ind_new) < 50
                continue
            else %if the unique object is > 50pixels then make it a seperate mask and update nIC
                cat(3, mask_cell_temp, zeros(size(mask_cell_temp(:,:,1))));
                nIC = nIC+1;
                if i == 1
                    %adds the new mask to the temp mask
                    mask_cell_temp(find(mask_cell(:,:,ic)== i),ic) = 1;
%                     ica_new(:,ic) = icasig(:,ic);
                else
                    %converts all non-1 values to ones and add the new mask to the end of the temp mask
                    mask_cell_temp(find(mask_cell(:,:,ic)== ind_rem(i)),nIC) = 1;
%                     ica_new(:,nIC) = icasig(:,ic);
                end
            end
        end
    else
        mask_cell_temp(find(mask_cell(:,:,ic)),ic) = 1;
%         ica_new(:,ic) = icasig(:,ic);
    end
end

%sort mask_cell according to mask size
[~, mask_cell_sizes_ind] = sort(squeeze(sum(sum(mask_cell_temp,2),1)), 'descend');
mask_cell_temp = mask_cell_temp(:,:,mask_cell_sizes_ind);

%finds overlapping pixels of ROIs and based on amount of overlap decides whether
%to group them or to split them- if splitting, then overlapping pixels are
%eliminated from both ROIs
nIC = size(mask_cell_temp,2);
mask_overlap = zeros(1,npw*nph);
mask_all = zeros(1,npw*nph);
min_size = 50;

for ic = 1:nIC
    ind_new = find(mask_cell_temp(:,ic))';
    %if there are more than 50 pixels in this mask
    if length(ind_new) >= min_size
        ind_old = find(mask_all); %finds all pixels (vectorized) serially updated full mask
        overlap = ismember(ind_old,ind_new);  %labels pixels of old mask which overlap with the current roi
        ind_both = find(overlap); %indeces of the overlapping pixels in terms of ind_old
        %if there are overlapping pixels..
        if length(ind_both)>0
%             write in portion to add a mask if another mask is completely ecnompassed by the new mask
%             if length(ind_new) - length(mask_all(ind_old(ind_both))) > min_size
%                 mask_all(ind_old(ind_both)) = 0; %converts all the overlapping pixels to 0
%                 mask_all(ind_new) = ic;
%             end

           %find the unique values for the overlapping pixels. ie. overlapping groups belonging to seperate masks
           ic_match = unique(mask_all(ind_old(ind_both)))';  
           % find how many overlapping pixels belong to each old mask and those maskIDs
           ic_match(:,2) = 0;
           for overlap_mask = 1:size(ic_match,1)
               ic_match(overlap_mask,2) = length(find(mask_all(ind_old(ind_both))==ic_match(overlap_mask,1))); 
           end
           %sort ic_match according to the amount of overlapping pixels
           [~,ic_match_sort_ind] = sort(ic_match(:,2), 'descend');
           ic_match = ic_match(ic_match_sort_ind,:);
           
           %determine if the mask with the most overlap mostly encompasses the new mask. If so
           %then label the new mask with that ID and remove portions of overlap with other masks.
           if ic_match(1,2)/length(ind_new) > 0.6 
               %assign new mask the old mask's ID and remove regions which overlap with other masks
               ind_old2 = find(mask_all~=ic_match(1,2) & mask_all~=0); %indeces in terms of full frame of all masks except the one which had overlap
               overlap2 = ismember(ind_old2, ind_new);  %indeces in terms of ind_old2
               ind_both2 = ind_old2(find(overlap2)); %indeces of overlapping pixels (except for the most overlapping) in terms of full frame
               ind_new( find(ismember(ind_new, ind_both2)) ) = [];  % remove overlapping pixels from the new mask
               mask_all(ind_new) = ic_match(1,2);
           else  %otherwise just label the new mask and get rid of overlapping portions
                mask_all(ind_new) = ic; %min(mask_all(ind_old(ind_both))); %looks at all the pixels from the new mask and sets them equal to ic
                mask_all(ind_old(ind_both)) = 0; %converts all the overlapping pixels to 0
                mask_overlap(ind_old(ind_both)) = 1;
           end
        else %if there are 0 overlapping pixels between the current roi and the serially updated mask then add this roi to the mask
            mask_all(ind_new) = ic;
        end
    end
end

mask_all = reshape(mask_all, npw, nph);
%mask_all = bwlabel(mask_all);
mask_final = reshape(mask_all, 1,npw*nph);

% removes ICs smaller than 200 pixels, renumbers ROIs so in continuous ascending order
% start_ic = 1;
% mask_final = zeros(size(mask_all));
% mask_all2 = mask_all;
% for ic = 1:max(mask_all2,[],2)
%     ind = find(mask_all2==ic);
%     if length(ind)<200
% %         mask_overlap(mask_all==ic) = 1;
%         mask_all2(ind) = 0;
%         ic
%     end
%     ind = find(mask_all2==ic);
%     if ~isempty(ind)
%         mask_final(ind)=start_ic;
%         start_ic= start_ic+1;
%     end
% end
% mask_final = reshape(mask_all2, 1,npw*nph);


end