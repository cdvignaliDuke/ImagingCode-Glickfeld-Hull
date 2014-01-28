outDir = fullfile(base,mouse,date,'analysis');
seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
load(fullfile(outDir,seqfile));

stack_sorted = zeros(size(stack) ,'uint16');

start = 1;
for iCond = 1:nCond;
    nRep = length(Big_Seqposition(iCond).ind);
    for iRep = 1:nRep;
        ind = Big_Seqposition(iCond).ind(iRep);
        stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)) = stack(:,:,1+((ind-1)*((nOFF+nON)/nPlanes)):((nOFF+nON)/nPlanes)+((ind-1)*((nOFF+nON)/nPlanes)));
        start = start+((nOFF+nON)/nPlanes);
    end
end

nblanks = length(Big_Seqposition(end).ind);
for iblank = 1:nblanks;
    ind = Big_Seqposition(end).ind(iblank);
    stack_sorted(:,:,start:start-1+((nOFF+nON)/nPlanes)) = stack(:,:,1+(ind-1)*((nOFF+nON)/nPlanes):((nOFF+nON)/nPlanes)+(ind-1)*((nOFF+nON)/nPlanes));
    start = start+((nOFF+nON)/nPlanes);
end
clear('stack');

stim_reps = zeros(1,nCond+1);
for iCond = 1:nCond+1;
    stim_reps(1,iCond) = length(Big_Seqposition(iCond).ind);
end

fn = fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_sorted.tif']);
writetiff(stack_sorted,fn);

fn2 = fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_reps.mat']);
save(fn2, 'stim_reps');

