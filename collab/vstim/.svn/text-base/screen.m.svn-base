function s = screen(monitor,varargin)
%SCREEN Creates screen structure
% S = SCREEN(MONITOR,VARARGIN)
%
%

if nargin < 1 
    monitor = 'vx922';
end

s = struct;

% here add your monitor type
switch monitor
    case 'vx922'
        s.Type = 'Viewsonic VX922';
        s.Size = [30.0 38.0]; % cm
        s.NominalResolution = [480 640];
        s.RefreshRate = 60; % Hz;
    case 'pixels'
        s.Size = [10 10]; % cm
        s.NominalResolution = [64 64];
        s.RefreshRate = 60; % Hz;
    otherwise
        error('unknown monitor type');
end

%% default arguments
defaultopts = {'UnderSamplingFactor',[2],'Distance',[15],'DecimationRatio',[1]};

for iarg = 1:2:length(defaultopts)
    s.(defaultopts{iarg}) =defaultopts{iarg+1};
end

for iarg = 1:2:length(varargin)
    s.(varargin{iarg}) =varargin{iarg+1};
end
   
s.Resolution = s.NominalResolution ./ s.UnderSamplingFactor ;
s.PixelSize = s.Size./s.Resolution; % cm

s.xpix = [0:s.Resolution(2)-1]-s.Resolution(2)/2;
s.ypix = [0:s.Resolution(1)-1]-s.Resolution(1)/2;

s.ydeg = pix2deg(s.ypix,s,1);
s.xdeg = pix2deg(s.xpix,s,2);
 
[s.X,s.Y]=meshgrid(s.xdeg,s.ydeg);

return;
