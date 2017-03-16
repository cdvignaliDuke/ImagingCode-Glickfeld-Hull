function mask_buffer = imCellBuffer(cellmask,buffersize)
% creates a buffer mask (FOV) around each selected cell 
%       cellmask = FOV of cell ROIs
%       buffersize = width of ring around each cell

seBR = strel('disk', buffersize); %buffer ring

siz = size(cellmask);
buffer = zeros(siz(1),siz(2),length(unique(cellmask))-1);
for i = 1:(length(unique(cellmask))-1)
    tempbuffering = cellmask == i;
    tempNPcells = cellmask;
    tempNPcells(tempNPcells == i) = 0; %dilate other cells so that you can subtract overlap
    cellBuffer = imdilate(tempNPcells,seBR);
    tempmask = imdilate(tempbuffering,seBR)-tempbuffering;
    tempmask = tempmask & ~cellBuffer;
    buffer(:,:,i) = tempmask;
end

cellsMaskCell = unique(cellmask);
B = ones(size(buffer));
for i = 1:(length(unique(cellmask))-1)
    B(:,:,i) = buffer(:,:,i)*i;
end
cellsWithBuffer = unique(B);

missingcells = setdiff(cellsMaskCell,[0; cellsWithBuffer]);

if ~isempty(missingcells)
    tempmask = zeros(size(cellmask));
    for i = missingcells
        buffer(:,:,i) = max(buffer(:,:,cellsWithBuffer),[],3);
    end
    disp('Cells are missing buffer, added all other buffer')
end

mask_buffer = buffer;

end


