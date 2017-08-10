function data = LoadSBX_imgOnly(mouse,date,ImgFolder,fName,varargin)


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
if ~isempty(varargin)
    nframes = varargin{1};    
else
    nframes = info.config.frames;
end

data = sbxread(fName,0,nframes);
if size(data,1) > 1
    data = data(1,:,:,:);
end
data = squeeze(data);
end