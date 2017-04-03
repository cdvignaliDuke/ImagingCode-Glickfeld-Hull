function [resp_match nTrials] = catIndexedTrialsCell(resp_cell,ind_cell)

nTrials = sum(cellfun(@(x) size(x,2),ind_cell));

if length(ind_cell) > 1
        resp_match_cell = cellfun(@(x,y) x(:,:,y), resp_cell, ind_cell,'unif', 0);

        resp_match = cat(3,resp_match_cell{1},resp_match_cell{2});
    
else    
    ind = ind_cell{1};
    ind = ind{1};
    if isempty(ind)
        resp_match = [];
    else
        r = resp_cell{1};
        
        resp_match = r(:,:,ind);
    end
end
end