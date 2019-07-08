function pix = stackGetPixVals(stack,mask)
    nc = max(mask(:));
    [nRows nCols nFrames]=size(stack);
    stackCols = reshape(stack, [nRows*nCols, nFrames]);
    labelCols = reshape(mask, [nRows*nCols, 1]);

    pix = cell(1,nc);
    for ic = 1:nc
        ind = labelCols == ic;
        pix{ic} = double(stackCols(ind,:));
    end
end