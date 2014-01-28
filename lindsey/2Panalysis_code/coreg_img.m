function [img_coreg] = coreg_img(img,target);
% The goal of this script is to take in an avi, and register that to a
% previously identified source, output an array containing this mov and
% save it to tiff

%% When clicking on cells to match up target and your given run, best
% results occur when you click on the cells that your cellsort has already
% identified

clear('input_points');
clear('base_points');

img(find(img>1e5)) = 0;
img = (img./max(max(abs(img)))); % files can't have any depth

target(find(target>1e5)) = 0;
target = (target./max(max(abs(target)))); 

[input_points, base_points] = cpselect(double(img),double(target),'Wait', true) 

%%
sz_target  = size(target);
mytform    = maketform('affine',input_points(1:3,:), base_points(1:3,:));
registered = imtransform(double(img),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);

 
%% Actual registration

sz = size(img);
img_coreg = [];
img_coreg(sz(1),sz(2)) = 0;
img_coreg = imtransform(img,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);