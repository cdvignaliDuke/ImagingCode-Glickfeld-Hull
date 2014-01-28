fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_reps.mat']);
load(fn);

stack_pn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_sorted.tif']);
stack = readtiff(stack_pn);

siz = size(stack);
stack_avg = zeros(siz(1), siz(2), nON+nOFF, nCond+1);

start = 0;
for iCond = 1:nCond+1;
    for iFrame = 1:(nOFF+nON)
    stack_avg(:,:,iFrame,iCond) = mean(stack(:,:,iFrame+start:nOFF+nON:start+((nOFF+nON)*stim_reps(iCond))),3);
    end
    start = (nOFF+nON)*stim_reps(iCond)+start;
end
stack_avg = uint16(stack_avg);

stack_uncorr = zeros(siz,'uint16');
start = 1;
for iCond = 1:nCond;
    nRep = stim_reps(iCond);
    for it = 1:nRep
        for iframe = 1:nOFF+nON;
            stack_uncorr(:,:,start) = stack(:,:,start)- stack_avg(:,:,iframe,iCond);
            start = start+1;
        end
    end
end

blank_avg = zeros(240, 256, 1);
blank_avg = mean(stack(:,:,start:end),3);
blank_avg = uint16(blank_avg);

nblanks = stim_reps(end);
for it = 1:nblanks*(nOFF+nON);
    stack_uncorr(:,:,start) = stack(:,:,start)- blank_avg;
    start = start+1;
end

uncorr_av = squeeze(mean(stack_uncorr,3));
uncorr_df = bsxfun(@minus, double(stack_uncorr), uncorr_av);
uncorr_tc = squeeze(mean(mean(uncorr_df,2),1));

stack_uncorr_flash = uint16(uncorr_df);
ind = find(uncorr_tc>1.5);
stack_uncorr_flash(:,:,ind) = [];

fn_out = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_uncorr.tif']);
writetiff(stack_uncorr_flash, fn_out);
clear stack_uncorr
clear stack_uncorr_flash
clear uncorr_df

