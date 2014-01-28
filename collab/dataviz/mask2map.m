function map = mask2map(mask,values)
%MASK2MAP
%MAP = MASK2MAP(MASK,VALUES)

ncells = max(unique(mask(:)));

map = zeros(size(mask));

for iCell = 1:ncells
    sel = find(mask==iCell);
    map(sel) = values(iCell);
end

return;
