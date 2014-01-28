function stack_writestacktotif(stack, outName)
%
%$Id: stack_writestacktotif.m 97 2008-03-17 20:14:13Z histed $

[nRows,nCols,nFrames,nPlanes] = size(stack);

for iF = 1:nFrames
    tFrame = permute(stack(:,:,iF,:), [1 2 4 3]); % no scaling now
    if iF == 1
        writeStr = 'overwrite';
    else
        writeStr = 'append';
    end        
    imwrite(tFrame, outName, ...
            'WriteMode', writeStr, ...
            'Compression', 'none');

    if mod(iF, 20) == 0
        fprintf(1, '%s: wrote %5d of %5d frames\n', mfilename, iF, nFrames);
    end
end


