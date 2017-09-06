function [data,nframes] = loadsbx_choosepmt(pmt,mouse,expdate,imgfolder,fname,varargin)

%% Set current directory to imaging data location
fdir = fullfile('Z:\data\', mouse, 'two-photon imaging',expdate,imgfolder);
cd(fdir);

%% load sbx info file
imgMatFile = [fname '.mat'];
load(imgMatFile);

%% get number of frames
if ~isempty(varargin)
    nframes = varargin{1}; 
else
    nframes = info.config.frames;
end

%% load data
data = sbxread(fname,0,nframes);
data = squeeze(data(pmt,:,:,:));
data = squeeze(data);
end