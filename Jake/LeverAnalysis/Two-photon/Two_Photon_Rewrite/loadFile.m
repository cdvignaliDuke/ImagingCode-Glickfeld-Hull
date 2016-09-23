function [data, skip_run] = loadFile(i)
file_info
data_dir = fullfile('\\crash\data\home\jake\Data\2P_imaging',[date{i} '_' mouseID{i}], mouseID{i});

config_fn = dir(fullfile(data_dir,['*' runID{i} '.mat']));
img_fn   = dir(fullfile(data_dir,['*' runID{i} '.sbx']));
[~,img_fn,~] = fileparts(img_fn.name);
if size(img_fn,1) ~= 0
    skip_run = 0;
    cd(data_dir);
    load(config_fn.name);
    nframes = info.config.frames;
    disp('loading sbx image file');
    data = squeeze( sbxread(img_fn,0,nframes) );
else
    skip_run = 1;
    data = [];
end