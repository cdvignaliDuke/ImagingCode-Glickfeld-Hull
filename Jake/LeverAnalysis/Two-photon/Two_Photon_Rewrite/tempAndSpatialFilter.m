function [data, varargout] = tempAndSpatialFilter(data,fps,varargin)

% pnrf = rangefilt(data, true([3 3 5])); % 6.5 minutes if [3 3 5], 31 sec if [3 3 3]
% data = bsxfun(@minus, data, mean(pnrf,3,'native')); % 2min
if nargin < 2
	fps = 20;
	fnyq = fps/2;
end
if nargin < 3
	
   % FIR FILTER
	n = 50;
	fstop = 5; %Hz
	wstop = fstop/fnyq;
	
	% DESIGNED FILTER
	d = designfilt('lowpassfir','SampleRate',fps, 'PassbandFrequency',fstop-.5, ...
		'StopbandFrequency',fstop+.5,'PassbandRipple',0.5, ...
		'StopbandAttenuation',65,'DesignMethod','kaiserwin');%could also use butter,cheby1/2,equiripple
else
	d = varargin{1};
end
data = temporalFilter(data,d);

	function dmat = temporalFilter(dmat,d)
		[phi,~] = phasedelay(d,1:5,fps);
		phaseDelay = mean(phi(:));
		h = d.Coefficients;
		h = double(h);
		filtPad = ceil(phaseDelay*4);
		% APPLY TEMPORAL FILTER
		sz = size(dmat);
		npix = sz(1)*sz(2);
		nframes = sz(3);
		% sdata = fftfilt(  (h), double( reshape(  (data), [npix,nframes])' ));
		sdata = filter( h, 1, double( cat(3, flip(dmat(:,:,1:filtPad),3),dmat)), [], 3);
		dmat = uint16(sdata(:,:,filtPad+1:end));
	end

if nargout > 1
	varargout{1} = d;
end

end