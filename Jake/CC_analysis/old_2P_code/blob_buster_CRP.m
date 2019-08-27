function mask_cell_out = blob_buster_CRP(mask_cell)

% get dimension sizes and allocate memory
[npw, nph, nIC] = size(mask_cell);
mask_cell_out = zeros(npw, nph, nIC);

%look at individual ICs and separate blobs which are not touching
for ic = 1:nIC
    if length(unique(bwlabel(mask_cell(:,:,ic))))>2 %if there is more than one object in this mask..
        this_mask_bw = bwlabel((mask_cell(:,:,ic)));
        for blob_num = unique(bwlabel(mask_cell(:,:,ic)))'
            if blob_num==0
                continue
            end
            blob_ind = find(this_mask_bw==blob_num); %mask_cell(find(this_mask_bw==blob_num));
            temp_mask = zeros(npw,nph);
            temp_mask(blob_ind) = 1; 
            if blob_num == 1
                mask_cell_out(:,:,ic) = temp_mask; 
            elseif blob_num > 1
                mask_cell_out(:,:,end+1) = temp_mask; 
            end
        end
    else
         mask_cell_out(:,:,ic) = mask_cell(:,:,ic);
    end
end
