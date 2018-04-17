function [data, xc] = correctMotion(data)
sz = size(data);
nFrames = sz(3);
% if nargin < 2
% 	
% 	prealign.n = 0;
% end
gpuVid = gpuArray(data);

% if ~isfield(prealign, 'template')
% 	templateFrame = im2single(gpuVid(:,:,1));
% else
% 	templateFrame = gpuArray(prealign.template);
% end

N = nFrames;
xc.cmax = zeros(N,1);
xc.xoffset = zeros(N,1); 
xc.yoffset = zeros(N,1);

% ESTIMATE IMAGE DISPLACEMENT USING NORMXCORR2 (PHASE-CORRELATION)
for k = 2:N
	movingFrame = gpuVid(:,:,k);
    templateFrame = gpuVid(:,:, k-1);
	c = normxcorr2(templateFrame, movingFrame);
	
	% RESTRICT VALID PEAKS IN XCORR MATRIX
% 	if isempty(validMaxMask)
% 		validMaxMask = false(size(c));
% 		validMaxMask(offsetShift-maxOffset:offsetShift+maxOffset, offsetShift-maxOffset:offsetShift+maxOffset) = true;
% 	end
% 	c(~validMaxMask) = false;
% 	c(c<0) = false;
	
	% FIND PEAK IN CROSS CORRELATION
	[cmax, imax] = max(abs(c(:)));
	[ypeak, xpeak] = ind2sub(size(c),imax(1));
	xoffset = sz(2) - xpeak;
	yoffset = sz(1) - ypeak;
	dx = gather(xoffset);
	dy = gather(yoffset);
	% APPLY OFFSET TO TEMPLATE AND ADD TO VIDMEAN
	adjustedFrame = imtranslate(gather(movingFrame), [dx, dy]);
% 	nt = prealign.n / (prealign.n + 1);
% 	na = 1/(prealign.n + 1);
% 	templateFrame = templateFrame*nt + gpuArray(adjustedFrame*na);
% 	prealign.n = prealign.n + 1;
	xc.cmax(k) = gather(cmax);
	xc.xoffset(k) = dx;
	xc.yoffset(k) = dy;
	
	% APPLY OFFSET TO FRAME
	
	data(:,:,k) = adjustedFrame;
	
end

end