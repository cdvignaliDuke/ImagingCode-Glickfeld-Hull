function [sort_mat, sort_ind] = sortCells4Heatmap(mat,trans_win)

    if ~isempty(mat)
        resp = nanmean(mat(trans_win,:),1);
        [resp_sort, sort_ind] = sort(resp);
        sort_mat = mat(:,fliplr(sort_ind))';
    else    
        sort_mat = [];
        sort_ind = [];
    end

end