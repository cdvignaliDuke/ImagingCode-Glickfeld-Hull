function [data, skip_run, img_fn] = loadFile(i,rID)
file_info_CRP
% data_dir = fullfile('\\crash\data\home\jake\Data\2P_imaging',[date{i} '_' mouseID{i}], mouseID{i});
%data_dir = fullfile('Z:\home\jake\Data\2P_imaging',[date{i} '_' mouseID{i}], mouseID{i});
data_dir = fullfile('Z:\Data\2P_imaging',[date{i} '_' mouseID{i}], mouseID{i});
config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
img_fn   = dir(fullfile(data_dir,['*' runID{rID} '.sbx']));

if size(img_fn,1) ~= 0
    [~,img_fn,~] = fileparts(img_fn.name);
    skip_run = 0;
    cd(data_dir);
    load(config_fn.name);
    if strcmp('img043_000_000', img_fn)
        nframes = 74463;
    else
        nframes = info.config.frames;
    end
    fprintf('loading sbx image file: %s %s\n', date{i}, img_fn);
    data = squeeze( sbxread(img_fn,0,nframes) );
else
    skip_run = 1;
    data = [];
end