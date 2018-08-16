function [nframes, skip_run, img_fn] = get_nFrames_temp(data_dir, runID, startframe, nframes)
% file_info
% data_dir = fullfile('\\crash\data\home\jake\Data\2P_imaging',[date{i} '_' mouseID{i}], mouseID{i});

config_fn = dir(fullfile(data_dir,['*' runID '.mat']));
img_fn   = dir(fullfile(data_dir,['*' runID '.sbx']));

if size(img_fn,1) ~= 0
    [~,img_fn,~] = fileparts(img_fn.name);
    skip_run = 0;
    cd(data_dir);
    load(config_fn.name);
    if strcmp('img043_000_000', img_fn)
        nframes = 74463;
    elseif ~isempty(strfind(data_dir, '170921_img044'))
        nframes = 100000;
    elseif ~isempty(strfind(data_dir, '171227_img067'));
        nframes = 73000;
    elseif ~isempty(strfind(data_dir, '180104_img067'));
        nframes = 87000;
    elseif isempty(nframes)
        nframes = info.config.frames;
    end
    
    
else
    skip_run = 1;
    nframes = [];
end