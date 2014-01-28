%average all frames of stack
outDir = fullfile(base,mouse,date,'analysis');
stack_sorted = readtiff(fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_sorted.tif']));

siz = size(stack_sorted);
stack_avg = zeros(siz(1), siz(2), nON+nOFF, nCond+1);
stack_avg_base = zeros(siz(1), siz(2), nCond+1);
stack_avg_resp = zeros(siz(1), siz(2), nCond+1);
start = 0;
for iCond = 1:nCond+1;
    for iFrame = 1:(nOFF+nON)
    stack_avg(:,:,iFrame,iCond) = mean(stack_sorted(:,:,iFrame+start:nOFF+nON:start+((nOFF+nON)*stim_reps(iCond))),3);
    end
    stack_avg_base(:,:,iCond) = mean(stack_avg(:,:,pre_win(1):pre_win(2),iCond),3);
    stack_avg_resp(:,:,iCond) = mean(stack_avg(:,:,post_win(1):post_win(2),iCond),3);
    stack_dF_all(:,:,iCond) = (stack_avg_resp(:,:,iCond)-stack_avg_base(:,:,iCond))./stack_avg_base(:,:,iCond);
    start = (nOFF+nON)*stim_reps(iCond)+start;
end

clear('stack_avg');
clear('stack_sorted');
fn_out = fullfile(outDir, [date '_' mouse '_run' num2str(userun) '_stack_dF_all.tif']);
writetiff(stack_dF_all, fn_out);