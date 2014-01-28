function [mov_coreg] = coreg_avi(mov,coreg_target, newdir);
% The goal of this script is to take in an avi, and register that to a
% previously identified source, output an array containing this mov and
% save it to tiff

%% When clicking on cells to match up target and your given run, best
% results occur when you click on the cells that your cellsort has already
% identified

% Master folder from which you will look at reconstructed images for both
% target and files you want to coregister
startupMA
PWD_input = 'E:\twophoton_data\2photon\raw';

% Example file you want to coregister
% file = 'NJ3\130719_NJ3\130719_NJ3_run2';
file = newdir;
expt = frGetExpt(file);
file_name = strtok(expt.name(end:-1:1), filesep);
file_name = file_name(end:-1:1);
if ~exist(expt.dirs.regrootpn)
    mkdir(expt.dirs.regrootpn);
end
if ~exist(expt.dirs.analrootpn)
    mkdir(expt.dirs.analrootpn);
end

% File you want to coregister to
% file_target = 'NJ3\130708_NJ3\130708_NJ3_run5';
file_target = coreg_target;
expt_coreg = frGetExpt(coreg_target);

%Place you want to store new ori_dff_tiff_mean and the rotation coordinates
% PWD_output = 'E:\twophoton_analysis\NJ3\130719_NJ3\130719_NJ3_run2\';
PWD_output = expt.dirs.analrootpn;

% AVG = double(readtiff([PWD_input,file])); % used to have ,1
target = double(load_tiff_folder([PWD_input,file_target]));
target = max(target,[],3);
% target = mean(target,3);
% target = target(:,:,1);
% mov = double(load_tiff_folder([PWD_input,file]));
AVG = max(mov,[],3);
% AVG = mean(mov,3);
% AVG = mov(:,:,1);

target0 = target;
AVG0 = AVG;

size(AVG)
size(target)
sz14 = size(AVG);
clear('input_points');
clear('base_points');

AVG(find(AVG>1e5)) = 0;
% AVG = (AVG./max(max(abs(AVG))))./2+.3;
AVG = (AVG./max(max(abs(AVG)))); % AVG can't have any depth

target(find(target>1e5)) = 0;
% target = (target./max(max(abs(target))))./2+.3; % Used to be plus 0.5
target = (target./max(max(abs(target)))); % Used to be plus 0.5
% cpselect(double(AVG),double(target)) 

[input_points, base_points] = cpselect(double(AVG),double(target),'Wait', true) 

%%
sz_target  = size(target);
mytform    = maketform('affine',input_points(1:3,:), base_points(1:3,:));
registered = imtransform(double(AVG0),mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);

%% NEED TO INSERT PROPER NAMING FORMAT
file_coordinates = [PWD_output,'\', file_name, '_correg_coordinates','.mat']
save(file_coordinates,'input_points','base_points','target','registered','AVG','file','mytform');
 
%% Actual registration

sz = size(mov);
mov_coreg = [];
mov_coreg(sz(1),sz(2),sz(3)) = 0;
mov_coreg = imtransform(mov,mytform,'XData',[1 sz_target(2)],'YData',[1 sz_target(1)]);

% writetiff(mov_coreg,[expt.dirs.regrootpn,'\', file_name],'double');