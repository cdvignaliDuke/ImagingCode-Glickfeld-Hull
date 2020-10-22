% get average image from example 2P session
% script for generating useable image to import in coreldraw
% this is fucking stupid. I don't know why:
% if import jpeg image into croeldraw: resolution is lower 
% can't import tiff into coreldraw
% so for the first 2 figures in figure 2, I used copy figure in matlab
% figure. 
% if copy figure right after imageshow, figure in coreldraw is wierd
% so in here we need to plot a black circle/dot on the black background and
% then copy figure into coreldraw, this works 
% fucking stupid

% sessions = '190603_img1025'; 
% image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 
% image_analysis_dest = [image_analysis_base, sessions, '\' 'getTC\'];
% eg = figure; imshow([image_analysis_dest 'AVG_' sessions '_000_rgstr_tiff_0_60000_50_ref42_' 'jpeg.jpg']); hold on;
% plot(1,1,'.','color','k');
% % bound1 = [148 149 150 151 152; 605 606 607 608 609];
% % plot(bound1(:,2),bound1(:,1),'.','color','k'); 
% savefig('Z:\Analysis\figures\figure2_2Prun\eg_ref_190603_img1025_addblack');

%%
% sessions = '191115_img1041'; 
% %ID = '1016';
% image_analysis_base    = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\'; 
% image_analysis_dest = [image_analysis_base, sessions, '\' 'getTC\'];
% eg = figure; imshow([image_analysis_dest 'AVG_' sessions '_000_rgstr_tiff_1_29999_50_ref23_' 'jpeg.jpg']); hold on;
% plot(1,1,'.','color','k');
% % bound1 = [148 149 150 151 152; 605 606 607 608 609];
% % plot(bound1(:,2),bound1(:,1),'.','color','k'); 
% savefig('Z:\Analysis\figures\figure4_airpuff_selfpace\eg_ref_191115_img1041_addblack');

%%
clear;
sessions = '190603_img1025'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 
image_analysis_dest = [image_analysis_base, sessions,'\'];
mask3D = load([image_analysis_dest 'getTC\' sessions '_000_nPCA200_mu0.3_nIC200_thresh96_coor0.8_mask3D_final.mat']);
mask3D = mask3D.mask3D;
threshold = -4;
spk_deconv_output = load([image_analysis_dest sessions,'_spk_deconvolve_threshold' num2str(threshold) '.mat']);
badPCs = spk_deconv_output.badPCs;  
eg_cells = [62,51,50,44,43,42,38,31,21,20,18,17,16,15,13,8];
allcells = 1:size(mask3D,3);
goodPCs = setdiff(allcells,badPCs);
other_cells = setdiff(goodPCs,eg_cells);

eg = figure; imshow([image_analysis_dest 'getTC\' 'AVG_' sessions '_000_rgstr_tiff_0_60000_50_ref42_jpeg.jpg']); hold on;
for i  = 1:length(other_cells)
    bound = cell2mat(bwboundaries(mask3D(:,:,other_cells(i))));
    %randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'color',[0.6706    0.8510    0.9137],'LineWidth',0.2); hold on;
end


for i = 1:length(eg_cells)
    bound = cell2mat(bwboundaries(mask3D(:,:,eg_cells(i))));
    %randcolor = rand(1,4);
    plot(bound(:,2),bound(:,1),'color',[1.0000 1.0000 0.7020],'LineWidth',1); hold on;
    text(mean(bound(:,2)),mean(bound(:,1)), num2str(i), 'color', [1.0000 1.0000 0.7020], 'FontSize', 10);
end
savefig('Z:\Analysis\figures\figure2_2Prun\eg_mask_190603_img1025_bicolor');

