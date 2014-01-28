ic4 = imscale(squeeze(ica_filters(4,:,:)))>.3;
figure; imagesc(ic4)

se1 = strel('disk',1);
ic4_close = imclose(ic4,se1);
ic4_open = imopen(ic4_close,se1);
figure; imagesc(ic4_open)
axon_mask = bwlabel(ic4_open);
figure; imagesc(axon_mask)
timeCourses = stackGetTimeCourses(stack,axon_mask);
figure; tcoffsetPlot(timeCourses)

nRoi = size(timeCourses,2);
r = zeros(nRoi);
for iroi1 = 1:nRoi;
    for iroi2 = 1:nRoi;
        r(iroi1, iroi2) = triu2vec(corrcoef(timeCourses(:,iroi1),timeCourses(:,iroi2)));
    end
end
figure;
imagesc(r);
colormap hot;
colorbar;

[order,groups] = corrsort(r);

figure;
imagesc(r(order, order));
colormap hot;

axon_fig = zeros(size(axon_mask));
for imask = 1:size(timeCourses,2);
    new_mask = find(axon_mask == imask);
    axon_fig(new_mask) = order(imask);
end
figure;
imagesc(axon_fig)

subplot(2,1,1);imagesc(lowpass(timeCourses(:,order)',[0 10]));colorbar
    