function [resp_mat_val, resp_mat_inv] = catExptAndNorm2BL(resp_cell_val, resp_cell_inv, bl_win)

    empty_ind_val = find(cellfun(@(x) isempty(x), resp_cell_val));
    empty_ind_inv = find(cellfun(@(x) isempty(x), resp_cell_inv));
    
    empty_ind = ones(size(resp_cell_val));
    empty_ind(unique([empty_ind_val, empty_ind_inv])) = 0;
    empty_ind = logical(empty_ind);
    
    temp_mat = cell2mat(resp_cell_val(empty_ind));
    bl = nanmean(temp_mat(bl_win,:),1);
    resp_mat_val = bsxfun(@minus, temp_mat, bl); 
    
    temp_mat = cell2mat(resp_cell_inv(empty_ind));
    bl = nanmean(temp_mat(bl_win,:),1);
    resp_mat_inv = bsxfun(@minus, temp_mat, bl); 

end