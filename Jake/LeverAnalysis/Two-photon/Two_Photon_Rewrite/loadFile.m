function [data, skip_run, img_fn] = loadFile(data_dir, runID, startframe, nframes)
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
    elseif ~isempty(strfind(data_dir, '170506_img90'));
        nframes = 115857;
    elseif isempty(nframes)
        nframes = info.config.frames;
    end
    fprintf('loading sbx image file: %s %s\n', data_dir, img_fn);
    if isempty(startframe)
        data = squeeze( sbxread([data_dir,img_fn],0,nframes) );
    else
        data = squeeze( sbxread([data_dir,img_fn],startframe,nframes) );
        
    end
    if length(size(data)) > 3
        data = squeeze(data(1,:,:,:));
    end
else
    skip_run = 1;
    data = [];
end