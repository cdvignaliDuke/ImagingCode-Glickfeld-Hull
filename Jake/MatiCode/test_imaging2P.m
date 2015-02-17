clear;

WRITE_VEDIO = 0;
BIN_SIZE = 1;   % if the ROI is large use large bin size (10)
PLOT_AIRPUFF = 1;
PLOT_RUN = 1;
% --- pcp2 mouse
 % image_dest  = 'Z:\Data\2P_imaging\150124img16\img16_000_000_cropped990.ome.tif';
 behave_dest =  'Z:\Data\2P_imaging\150126_img16\data-img16-150126-1802';
 % first_frame = 25; % In this day the first 25 frames were off
  first_frame = 1; % In this day the first 25 frames were off
  last_frame = 4950;
  Sampeling_rate = 15.5;

%  image_dest  = '../data/141126_img15/141126_img15_2_1_MMStack.ome.tif';
%  behave_dest =  '../data/141126_img15/data-141126-1604.mat';

%   image_dest  = '../data/141119_img9/141119_img9_1_MMStack.ome.tif';
%   behave_dest =  '../data/141119_img9/data-141119-1430.mat';
%   first_frame = 1; % In this day the first 25 frames were off
%   Sampeling_rate = 4;

% image_dest  = '../data/150109_img17/150109_img17_1.tif';
% behave_dest =  '../data/150109_img17/150109_img17_1.mat';
% first_frame = 1;
% last_frame = 1200;
% Sampeling_rate = 4;
% image_dest  = '../data/150109_img18/150109_img18_2.tif';
% behave_dest =  '../data/150109_img18/150109_img18_2.mat';
% 
% 
% first_frame = 1;
% last_frame = 200;
% Sampeling_rate = 4;

% image_dest  = '/Volumes/Promise RAID/mati/imaging_data/150109_img16/150109_img16_1.tif';
% behave_dest =  '/Volumes/Promise RAID/mati/imaging_data/150109_img16/150109_img16_1.mat';

%  image_dest  = '../data/150109_img16/150109_img16_1.tif';
%  behave_dest =  '../data/150109_img16/150109_img16_1.mat';
% 
% first_frame = 1;
% last_frame = 1200;
% Sampeling_rate = 4;

% image_dest  = '../data/150109_img18/150109_img18_3.tif';
% behave_dest =  '../data/150109_img18/150109_img18_3.mat';
% first_frame = 1;
% last_frame = 1200;
% Sampeling_rate = 4;



b_data = load(behave_dest);
[mouse_run, airpuff] = parse_behavior_for_running(b_data.input); % get running and airpuff times

info = imfinfo(image_dest);
if(~exist('first_frame', 'var'))
    first_frame = 1;
end
if(~exist('last_frame', 'var'))
    last_frame =  numel(info);
end

[ROI_x, ROI_y] = get_movie_ROI(info, image_dest, first_frame, last_frame); % get the ROI -  must be a rectangle

[img, sz] = get_movie_by_ROI(image_dest, info, ROI_x, ROI_y, BIN_SIZE, first_frame, last_frame);
if(WRITE_VEDIO)
    params.dest = image_dest(1:end-4);
    params.sampeling_rate =Sampeling_rate;
    params.sz = sz;
    params.speed = 1;
    params.remove_avg = 1;
    params.zscroe = 0;
    params.colormap  = 'jet';
    write_movie(img, params);
end
% calculate average and SD
avg = mean(img,2);

diviation = bsxfun(@minus, img, avg);
sq_d =    diviation.*diviation;
sd = sqrt(mean(sq_d,2));
avg_across_pix =  nanmean(diviation);

% --------- plot SD and average
figure; subplot(2,2,1); pcolor(reshape(avg,sz)); axis ij; shading flat;
title('mean');colorbar;
subplot(2,2,2); pcolor(reshape(sd,sz)); axis ij; shading flat;
title('SD');colorbar;

% calculate the covariance of the max SD point
[ val max_sd_inx] = max(sd);

max_div = diviation(max_sd_inx,:);
mul_with_max = bsxfun(@times, diviation, max_div);
cov_vec = mean(mul_with_max,2);
max_vals = img(max_sd_inx,:);
cc_vec = corr(max_vals', img');

subplot(2,2,3); pcolor(reshape(cov_vec,sz)); axis ij; shading flat;
title('COV with max SD');
colorbar;

subplot(2,2,4); pcolor(reshape(cc_vec,sz)); axis ij; shading flat;
title('CC with max SD'); colorbar;

figure;
plot(first_frame:last_frame, zscore(max_vals) ,'g'); hold on;
plot(first_frame:last_frame, zscore(avg_across_pix), 'k'); hold on;

if(PLOT_RUN)
    p_run  = mouse_run(1001:4000);
    p_run = p_run -nanmean(p_run);
    p_run = p_run./nanstd(p_run);
    plot(first_frame:3000, p_run,'r'); hold on;
end

if(PLOT_AIRPUFF)
    air_puff_time = find(airpuff >0);
    p_air_puff_time = air_puff_time(air_puff_time >= first_frame & air_puff_time <=last_frame);
    plot(p_air_puff_time, 2, '.g', 'markersize' ,10)
end


%  --- plot difference between  running and not running images
if(PLOT_RUN)
    % ------- plot trials with and without running - beware  of NaNs
    is_running_inx = find(mouse_run >= 1) -first_frame;
    is_running_inx = is_running_inx(is_running_inx >0 & is_running_inx <=last_frame-first_frame);
    not_running_inx = find(mouse_run == 0) -first_frame;
    not_running_inx = not_running_inx(not_running_inx> 0 & not_running_inx <= last_frame-first_frame);
    
    is_running = img(:,is_running_inx);
    not_running = img(:,not_running_inx);
    avg_run = mean(is_running, 2);
    avg_not_run = mean(not_running, 2);
    
    figure;
    max_v =  max([avg_run;avg_not_run]);
    subplot(2,2,1); pcolor(reshape(avg_run, sz));  title('running'); axis ij; shading flat; caxis([0 max_v]); colorbar;
    subplot(2,2,2); pcolor(reshape(avg_not_run, sz)); title('not running'); axis ij; shading flat;caxis([0 max_v]); colorbar;
    subplot(2,2,3); pcolor(reshape(avg_run-avg_not_run, sz)); title('run - not running'); axis ij; shading flat; colorbar;
    corr_with_run = corr(mouse_run(first_frame:last_frame), img', 'rows' , 'pairwise');
    subplot(2,2,4); pcolor(reshape(corr_with_run, sz)); title('correlation with run'); axis ij; shading flat; colorbar;
end


%  --- plot difference between  images with/without  airpuff.
if(PLOT_AIRPUFF)
    % -- plot trials with and without airpuff
    is_airpuff_inx = find(airpuff >= 1) -first_frame;
    is_airpuff_inx = is_airpuff_inx(is_airpuff_inx >0 & is_airpuff_inx <= last_frame-first_frame);
    not_airpuf_inx = find(airpuff == 0) -first_frame;
    not_airpuf_inx = not_airpuf_inx(not_airpuf_inx> 0 & not_airpuf_inx <= last_frame-first_frame) ;
    
    is_airpuff = img(:,is_airpuff_inx);
    not_airpuff = img(:,not_airpuf_inx);
    
    avg_air = mean(is_airpuff, 2);
    avg_not_air = mean(not_airpuff, 2);
    
    max_v =  max([avg_air;avg_not_air]);
    figure;
    subplot(2,2,1); pcolor(reshape(avg_air, sz));  title('airpuff'); axis ij; shading flat; caxis([0 max_v]); colorbar;
    subplot(2,2,2); pcolor(reshape(avg_not_air, sz)); title('no airpuff'); axis ij; shading flat;caxis([0 max_v]); colorbar;
    subplot(2,2,3); pcolor(reshape(avg_air-avg_not_air, sz)); title('airpuff - no airpuff'); axis ij; shading flat; colorbar;
    
end
