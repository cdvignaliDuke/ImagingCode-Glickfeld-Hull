fn = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn);

nReps = sum(stim_reps(1,1:nCond));
nblanks = stim_reps(1,end);

resp = zeros(nReps/nCond, nON+nOFF, nCond);
blanks = zeros(nblanks, nON+nOFF);

start = 1;
for icon = 1:nCond;
    for iRep = 1:stim_reps(1,icon);
    resp(iRep, :, icon) = roi_avg(1, start:start-1+nON+nOFF);
    start = start+nON+nOFF;
    end
end

for iblank = 1:nblanks;
    blanks(iblank, :) = roi_avg(1, start:start-1+(nON+nOFF));
    start = start+nON+nOFF;
end
resp_norm = bsxfun(@minus, resp,mean(resp(:,Pre_win(1):Pre_win(2),:),2));
resp_dF = 
blanks_norm = bsxfun(@minus, blanks,mean(blanks(:,Pre_win(1):Pre_win(2)),2));

fn = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_resp_norm.mat']);
save(fn, 'resp_norm', 'blanks_norm');