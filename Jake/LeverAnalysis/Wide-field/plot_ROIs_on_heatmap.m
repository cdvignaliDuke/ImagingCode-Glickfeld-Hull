function plot_ROIs_on_heatmap(avg_img, sz, b_data, day, cluster);
%plot the heatmap and overlay outlines of ROIs with numbers

%different sized figures for lever/no-lever and Cue-Reward pairing
if b_data.doFakeMouseSuccessOnly == 0;
    subplot(2,3,1); 
elseif b_data.doFakeMouseSuccessOnly == 1;
    subplot(1,3,1); 
end
imagesc(reshape(avg_img, sz(1), sz(2)));

%give an appropriate title based on exp type
if b_data.doFakeMouseSuccessOnly;
    title(['Cue-Reward Pairing: ', day]);
elseif b_data.input.doLever == 0;
    title(['CONTROL: no lever ' day]);
else
    title(day);
end

%apply shading and plot ROI outlines
shading flat; hold on;
for i=1:cluster.num_cluster
    line(  cluster.roi_position{i}(:,1),   cluster.roi_position{i}(:,2) ,'color', 'k', 'linewidth', 2)
    text(mean(cluster.roi_position{i}(:,1)),mean(cluster.roi_position{i}(:,2)), ...
        [num2str(i)], 'color', 'k', 'FontSize', 30);
end

return