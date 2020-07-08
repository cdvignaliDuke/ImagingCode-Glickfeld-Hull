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
sessions = '190603_img1025'; 
%ID = '1016';
image_source_base  = 'Z:\Data\2photon\'; %location of permanently stored image files for retreiving meta data
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; 
image_analysis_dest = [image_analysis_base, sessions, '\' 'getTC\'];
eg = figure; imshow([image_analysis_dest 'AVG_' sessions '_000_rgstr_tiff_0_60000_50_ref42_' 'jpeg.jpg']); hold on;
plot(1,1,'.','color','k');
% bound1 = [148 149 150 151 152; 605 606 607 608 609];
% plot(bound1(:,2),bound1(:,1),'.','color','k'); 
savefig('Z:\Analysis\figures\figure2_2Prun\eg_ref_190603_img1025_addblack');
