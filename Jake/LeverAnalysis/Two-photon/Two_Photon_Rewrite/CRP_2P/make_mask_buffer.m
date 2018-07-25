function mask_3D_buffed = make_mask_buffer(mask_3D);
%function for creating a buffer region around each mask. If there exist
%pixels from another mask within the buffer region of a current mask then
%both pixels are eliminated. 

mask_3D_buffed = mask_3D;
%get number of masks
nmasks = size(mask_3D,3);

for this_mask = 1:nmasks
    %get 2D mask for all masks except the current mask
    curr_mask_2D = sum(mask_3D(:,:,[[1:this_mask-1], [this_mask+1:end]]),3);
    curr_mask_2D = logical(curr_mask_2D);
    
    %find all pixels in the buffer zone of current mask (pixel group A)
    [cur_row, cur_col] = find(mask_3D(:,:,this_mask));
    
    %find all pixels=1 in 2D mask within the buffer zone of the current mask.
    %Define pixel group B as all pixels where other masks overlap with pixel group A.
    buffer_zone = mask_3D(:,:,this_mask);
    buffer_zone( sub2ind(size(buffer_zone), [cur_row+1], [cur_col+1]) ) = 1;
    buffer_zone( sub2ind(size(buffer_zone), [cur_row-1], [cur_col-1]) ) = 1;
    buffer_zone( sub2ind(size(buffer_zone), [cur_row+1], [cur_col-1]) ) = 1;
    buffer_zone( sub2ind(size(buffer_zone), [cur_row-1], [cur_col+1]) ) = 1;
    buffer_zone( sub2ind(size(buffer_zone), [cur_row], [cur_col+1]) ) = 1;
    buffer_zone( sub2ind(size(buffer_zone), [cur_row], [cur_col-1]) ) = 1;
    buffer_zone( sub2ind(size(buffer_zone), [cur_row+1], [cur_col]) ) = 1;
    buffer_zone( sub2ind(size(buffer_zone), [cur_row-1], [cur_col]) ) = 1;
    
    %store pixel indeces of pixels belonging to other masks within the current
    %masks buffer zone. Label as pixel group B.
    [B_row, B_col] = find((buffer_zone + curr_mask_2D) > 1 );
    
    %in final 3D mask make all those pixels = 0 (pixel group B). (exclude current mask)
    for mask_num = [[1:this_mask-1], [this_mask+1:nmasks]]
        lin_ind = sub2ind(size(mask_3D_buffed), [B_row], [B_col], repmat(mask_num, 1, length(B_col))');
        mask_3D_buffed(lin_ind) = 0;
    end
    
    %find all pixels within buffer zone of Pixel group B (pixel group C). Set all pixels on
    %the current mask within this buffer zone = 0 on the final 3D mask.
    group_C = zeros(size(curr_mask_2D));
    group_C( sub2ind( size(group_C), [B_row], [B_col]) ) = 1;
    
    group_C(sub2ind( size(group_C), [B_row+1], [B_col+1])) = 1;
    group_C(sub2ind( size(group_C), [B_row-1], [B_col-1])) = 1;
    group_C(sub2ind( size(group_C), [B_row-1], [B_col+1])) = 1;
    group_C(sub2ind( size(group_C), [B_row+1], [B_col-1])) = 1;
    group_C(sub2ind( size(group_C), [B_row], [B_col+1])) = 1;
    group_C(sub2ind( size(group_C), [B_row], [B_col-1])) = 1;
    group_C(sub2ind( size(group_C), [B_row-1], [B_col])) = 1;
    group_C(sub2ind( size(group_C), [B_row+1], [B_col])) = 1;
    
    [C_row, C_col] = find(group_C);
    
    lin_ind = sub2ind(size(mask_3D_buffed), [C_row], [C_col], repmat(this_mask, 1, length(C_col))'); 
    mask_3D_buffed(lin_ind) = 0;
    
    
end
return









