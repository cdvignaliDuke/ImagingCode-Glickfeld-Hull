function [img_mat_file, laser_power_vec] = get_laser_power_data(sub, rID);
% function for detecting if a 2P dataset used the laser toggle to turn off
% the laser during the ITI. Allows for motion registration and ICA/PCA to
% be limited to frames with laser power

%uses two files saved alongside the .sbx file. The info.frame variable in
%the .mat file logs the frame in which the laser turned on and off. The
%ttl_log variable in the _realtime file logs a zero for when the laser is
%off and a 1 for when it is on for each frame
file_info_CRP;
data_dir = fullfile('Z:\Data\2P_imaging',[date{sub} '_' mouseID{sub}], mouseID{sub});
config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
laser_power_fn = dir(fullfile(data_dir,['*' runID{rID} '_realtime.mat']));

%load imaging .mat file
if size(config_fn,1) ~= 0
    cd(data_dir);
    img_mat_file = load(config_fn.name);
else 
    img_mat_file = [];
end

%load realtime file 
if size(laser_power_fn,1) ~= 0
    load(laser_power_fn.name);
    laser_power_vec = ttl_log; 
else 
    laser_power_vec= [];
end

end







        
        
        
        