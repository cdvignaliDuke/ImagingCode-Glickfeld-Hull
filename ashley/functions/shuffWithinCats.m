function data_shuff = shuffWithinCats(data,categoryIDs,varargin)
% DATA should be organized as trials x cells, where trials will be the
% dimension that is randomized
% CATEGORYIDS should be a cell array where each cell is of equal length,
% identifing the binary category ID category a given trial belongs to
% varargin can be left empty for the possibility of seeding the
% randomization

    if ~isempty(varargin)
        rng(varargin{1})
    end
    
    [~,nc]=size(data);
    ncats = length(unique(categoryIDs));
    if min(categoryIDs) == 0
        categoryIDs = categoryIDs+1;
    end

    data_cats = cell(1,ncats);
    for icat = 1:ncats
        ind = categoryIDs == icat;
        data_cats{icat} = data(ind,:);
    end

    data_shuff = cell(size(data_cats));
    for icat = 1:ncats
        nt = size(data_cats{icat},1);
        d = data_cats{icat};
        for icell = 1:nc
            data_shuff{icat}(:,icell) = d(randperm(nt),icell);
        end
    end
end