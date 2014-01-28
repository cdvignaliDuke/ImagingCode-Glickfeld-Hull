fn = fullfile(outDir,[date '_' mouse '_run' num2str(userun) '_Big_Seqposition']);
load(fn);

figure;
subs = ceil(sqrt(size(roi_avg,2)));
col{1} = '-r';
col{2} = '-m';
col{3} = '-y';
col{4} = '-g';
col{5} = '-c';
col{6} = '-b';
for iRoi = 1:size(roi_avg,2);
    subplot(subs, subs, iRoi);
    start = 1;
    for nvar = 1:length(uvar);
        errorbar(uvar(:,1), resp_avg(iRoi,start:start+(length(uvar)-1)),resp_std(iRoi,start:start+(length(uvar)-1)),col{nvar});
        hold on;
        start = start+length(uvar);
    end
    hold on
    errorbar(-1,resp_avg(iRoi,end),resp_std(iRoi,end),'*k')
    xlim([-2 18])
end


for iRoi = 1:size(resp_avg,1);
    resp_peak_cond(iRoi) = max(resp_avg(iRoi,1:nCond),[],2);
    resp_avg_norm(iRoi,:) = resp_avg(iRoi,1:nCond)./max(resp_avg(iRoi,1:nCond),[],2);
end
resp_avg_norm(find(resp_avg_norm<0))=0;

resp_avg_norm_square = zeros(length(uvar), length(uvar), size(resp_avg,1));
for iRoi = 1:size(resp_avg,1);
    resp_avg_norm_square(:,:,iRoi) = reshape(resp_avg_norm(iRoi, :),[length(uvar) length(uvar)])';
end

figure;
subs = ceil(sqrt(size(resp_avg,1)));
for iRoi = 1:size(roi_avg,2);
    subplot(subs, subs, iRoi);
    imagesc(squeeze(resp_avg_norm_square(:,:,iRoi)));
    title(num2str(resp_peak_cond(iRoi)));
    xlabel('TF')
    ylabel('SF')
    colormap gray;
end


