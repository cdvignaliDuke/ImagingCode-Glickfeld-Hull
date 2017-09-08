[x,y] = size(clusters);

h = impoly;
junk_mask = (createMask(h));

mask = junk_mask & clusters;
mask = reshape(mask, x, y);
clusters_idx = reshape(clusters, x, y);
clusters_idx(mask) = 0;
clusters = reshape(clusters_idx, x, y);
figure;imagesc(clusters)

fig = figure; imagesc(clusters);
saveas(fig, 'clusters_wholeMovie_clean.fig');
save('clusters_wholeMovie_clean.mat', 'clusters');

bw = imerode(clusters, true(2));
se = strel('disk',1);
bw = imdilate(bw,se);

% figure
% for i = size(cSizeNorm,1)
%     
%     bar(1:length(cSizeNorm{i}, cSizeNorm{i}));
%     bar(length(cSizeNorm{i}+2) : , cSizeNorm{i}));
% end