fn = fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_Big_Seqposition']);
load(fn);

for iS = 1:length(Big_Seqposition)
    vars(iS)=Big_Seqposition(iS).TFSFetc(var);
end
uvars = unique(vars);
nvars = length(uvars);

figure;
subs = ceil(sqrt(size(roi_avg,2)));
for iRoi = 1:size(roi_avg,2);
    subplot(subs, subs, iRoi);
    errorbar(uvars, resp_avg(iRoi,1:nCond),resp_std(iRoi,1:nCond));
    hold on
    errorbar(-1,resp_avg(iRoi,end),resp_std(iRoi,end),'*k')
end
