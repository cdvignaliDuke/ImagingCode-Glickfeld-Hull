function data = loadsbx_pmt2(mouse,date,imgfolder,fName,varargin)

% Set current directory to imaging data location
fdir = fullfile('Z:\data\', mouse, 'two-photon imaging',expdate,imgfolder);
cd(fdir);
imgMatFile = [fName '.mat'];
load(imgMatFile);

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

data = sbxread(fName,0,nframes);
if size(data,1) > 1
    data = data(1,:,:,:);
end
% pmt = 2; %1 = green 2 = red
% data = squeeze(data(pmt,:,:,:,:));
data = squeeze(data);
end
end