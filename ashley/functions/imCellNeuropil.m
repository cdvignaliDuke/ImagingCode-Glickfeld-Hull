function f = imCellNeuropil(cellmask,buffersize,neuropilsize)
% creates a neuropil mask (FOV) around each selected cell 
%       cellmask == FOV of cell ROIs
%       buffersize = width of ring around each cell
%       neuropilsize = width of neuropil ring

seBR = strel('disk', buffersize); %buffer ring
seNPR = strel('disk', buffersize+neuropilsize); %neuropil ring

neuropil = zeros(size(cellmask));
for i = 1:(length(unique(cellmask))-1)
    tempmask = cellmask == i; %logical
    tempbufferring = cellmask == i;
    tempNPcells = cellmask;
    tempNPcells(tempNPcells == i) = 0; %dilate other cells so that you can subtract overlap
    cellBuffer = imdilate(tempNPcells,seNPR);
    tempmask = imdilate(tempmask,seNPR) - imdilate(tempbufferring,seBR);
    tempmask = tempmask & ~cellBuffer;
    neuropil(tempmask == 1) = i;
end

f = neuropil;

end


