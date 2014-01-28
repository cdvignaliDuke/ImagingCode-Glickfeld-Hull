function dec = tcDecimate(traces,r);
% DEC = TCDECIMATE(TRACES,R);

if ~isa(traces,'dougle')
    traces = double(traces);
end

[nSamples,nCells]=size(traces);

for iCell = 1:nCells;
   dec(:,iCell)= decimate(traces(:,iCell),r); 
end

return;