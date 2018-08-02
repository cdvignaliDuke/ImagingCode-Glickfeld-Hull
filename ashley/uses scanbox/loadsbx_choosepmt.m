function [data,nframes] = loadsbx_choosepmt(pmt,mouse,expdate,imgfolder,fname,varargin)

%% Set current directory to imaging data location
if length(varargin) > 1
%     rc = varargin{2};
%     fdir = fullfile(rc.ashleyData,mouse,'two-photon imaging',expdate,imgfolder);
    fdir = varargin{2};
else
    fdir = fullfile('\\crash.dhe.duke.edu\data\home\ashley\data', mouse, 'two-photon imaging',expdate,imgfolder);
end
cd(fdir);

%% load sbx info file
imgMatFile = [fname '.mat'];
load(imgMatFile);

%% get number of frames
if ~isempty(varargin)
    if ~isempty(varargin{1})
        nframes = varargin{1}; 
    else
        nframes = info.config.frames;
    end
else
    nframes = info.config.frames;
end

%% load data
data = sbxread(fname,0,nframes);
data = squeeze(data(pmt,:,:,:));
data = squeeze(data);
end