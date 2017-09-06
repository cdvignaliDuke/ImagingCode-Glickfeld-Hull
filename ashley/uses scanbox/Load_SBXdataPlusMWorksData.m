function [input, data, t] = Load_SBXdataPlusMWorksData(SubNum,date,time,mouse,ImgFolder,fName,varargin)
CD = 'Y:\home\andrew\Behavior\Data';
mworks = ['data-' 'i' SubNum '-' date '-' time]; 
load (fullfile(CD,mworks));

% Set current directory to imaging data location
CD = ['Z:\data\' mouse '\two-photon imaging\' date '\' ImgFolder];
cd(CD);
% CD = ['D:\Ashley_temp\' date '\' ImgFolder];
% cd(CD);
imgMatFile = [fName '.mat'];
load(imgMatFile);

% eyeMatFile = [fName '_eye.mat'];
% load(eyeMatFile);

%%
if nargout > 1
%new datasets
if ~isempty(varargin)
    nframes = varargin{1}; 
    if isempty(nframes)
        nframes = info.config.frames;
    end
else
    nframes = info.config.frames;
end

tic
data = sbxread(fName,0,nframes);
if size(data,1) > 1
    data = data(1,:,:,:);
end
t = toc;
% pmt = 1; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
data = squeeze(data);
end
end