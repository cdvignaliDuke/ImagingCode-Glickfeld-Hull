function outRaster = rasterflat(inRaster,trialdur)
% OUTRASTER = RASTERFLAT(INRASTER,TRIALDUR)

[nsamples,ncells]=    size(inRaster);

ntrials = nsamples / trialdur;

outRaster = zeros(trialdur,ncells*ntrials);

for iTrial = 1:ntrials
    rows = (iTrial-1)*trialdur+1:iTrial*trialdur;    
    cols = (iTrial-1)*ncells+1:iTrial*ncells;    
    outRaster(:,cols) =inRaster(rows,:);
end

return
