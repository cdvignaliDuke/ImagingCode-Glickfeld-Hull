function [piezo_data, frame_nums] = load_piezo(file_name)
%function for opening and reading the piezo data collected during 2P
%imaging sessions. 
%input: file_name should be a path and file name to a .ephys file
%outputs: piezo_data will be a timecourse of the piezo readings. 
%         frame_nums will be an index which tells you which imaging frame
%         the piezo data corresponds to. 

piezo_open = fopen(file_name);
piezo_data_mixed = fread(piezo_open, 'single');
piezo_data = piezo_data_mixed(2:2:end);
frame_nums = piezo_data_mixed(1:2:end);

return