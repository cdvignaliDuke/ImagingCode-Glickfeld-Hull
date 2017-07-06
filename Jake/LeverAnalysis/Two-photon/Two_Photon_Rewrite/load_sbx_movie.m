function [info, img, t] = load_sbx_movie(img_fn, data_dir, varargin)
% Set current directory to imaging data location
clear global
CD = [data_dir, '\', ];
cd(CD);
imgMatFile = [img_fn '.mat'];
load(imgMatFile);

%%
nframes = info.config.frames;
if ~isempty(varargin)
    num_good_frames = cell2mat(varargin);
    if num_good_frames > nframes
        num_good_frames=nframes;
    end
end

tic
img = sbxread(img_fn,0,num_good_frames);
t = toc;
img = squeeze(img);
info = 1;
end

