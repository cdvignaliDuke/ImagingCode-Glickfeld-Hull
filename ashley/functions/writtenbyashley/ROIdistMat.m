function distanceMatrix = ROIdistMat(mask_cell)
%% calculate the distance between ROIs in a cell mask, in pixels

[c ci] = unique(mask_cell(:));
nCells = length(c)-1;
[yCell xCell] = ind2sub(size(mask_cell),ci(2:end));

distanceMatrix = zeros(nCells);
c = zeros(1,nCells);
for iCell = 1:nCells
    for iC2 = 1:nCells
        if iC2 ~= iCell
            a = xCell(iCell) - xCell(iC2);
            b = yCell(iC2) - yCell(iCell);
            c(:,iC2) = sqrt((a^2)+(b^2));
        end
    end
    distanceMatrix(iCell,:) = c;
    c = zeros(1,nCells);
end


end