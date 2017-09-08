function roi = detectSingleFrameRois(data,info)
% INPUT:
%	Expects vid.cdata with cdata datatype = 'uint8'
% OUTPUT:
%	Array of objects of the 'RegionOfInterest' class, resembling a RegionProps structure
%	(Former Output)
%	Returns structure array, same size as vid, with fields
%			bwvid =
%				RegionProps: [12x1 struct]
%				bwMask: [1024x1024 logical]

% GET SAMPLE VIDEO STATISTICS AND DEFINE MIN/MAX ROI AREA
sz = size(data);
N = sz(3);
frameSize = sz(1:2);
% dsamp = getDataSample(data);
% stat.Min = min(data,[],3);
% stat.Std = std(double(data),1,3);
minRoiPixArea = 50; 
maxRoiPixArea = 500; 
maxRoiEccentricity = .92;%large value for elliptical shape
% maxPerimOverSqArea = 0.52; %  circle = 3.5449, square = 4 % Previousvalues: [6.5  ]
% minPerimOverSqArea = 3.0; % previously 3.5 PERIMETER / SQRT(AREA)
% INITIALIZE DYNAMIC SIGNAL THRESHOLD ARRAY: ~1 STD. DEVIATION OVER MINIMUM (OVER TIME)
% stdOverMin = 2; % formerly 1.2
% signalThreshold = gpuArray( stat.Min + uint16( stat.Std.*stdOverMin ));
% RUN A FEW FRAMES THROUGH HOTSPOT FINDING FUNCTION TO IMPROVE INITIAL SIGNAL THRESHOLD
% for k = fliplr(round(linspace(1,N,min(20,N))))
% 	[~, signalThreshold] = getAdaptiveHotspots(data(:,:,k), signalThreshold);
% end
% GET HOTSPOT-BASED ROIS (WHENEVER SIGNAL INTENSITY EXCEEDS PIXEL-SPECIFIC THRESHOLD)
% bwmask = false([frameSize N]);
if nargin<2
	info = [];
	frameNum = 1:N;
else
	frameNum = cat(1,info.frame);
end
frameROI = cell(N,1);
Ithreshold = 82; %85 for 0420_img92, 90 for 0427_img90, 
for k = 1:N
    k
% 	[bwmask(:,:,k), signalThreshold] = getAdaptiveHotspots(data(:,:,k), signalThreshold);
    bw = gpuArray(zeros(sz(1), sz(2)));
    diffImage = gpuArray(data(:,:,k));
%   bw = diffImage > signalThreshold;
		% changed from: bw = imclose(imopen( bw, S.disk6), S.disk4);
    bw (diffImage > 0.65*mean([max(prctile(diffImage,Ithreshold,1)) max(prctile(diffImage,Ithreshold,2))])) = 1;
%     bw = gather(bwmorph(bwmorph( bw, 'open'), 'majority'));
    bw = gather(bwmorph( bw, 'open'));
% 	bw = gather(bwmorph(bwmorph(bwmorph( bw, 'open'), 'shrink'), 'majority'));%4
%   bw = gather(bwmorph(bwmorph(bwmorph( bw, 'hbreak'), 'spur'), 'majority'));
    bwRP =  regionprops(bw,...
		'Centroid', 'Area', 'BoundingBox','PixelIdxList',...
		'Eccentricity', 'Orientation','Perimeter');
%     bwmask(:,:,k) = bw;
	bwRP = bwRP([bwRP.Area] >= minRoiPixArea);	%	Enforce MINIMUM SIZE
 	bwRP = bwRP([bwRP.Area] <= maxRoiPixArea);	%	Enforce MAXIMUM SIZE
	bwRP = bwRP([bwRP.Eccentricity] >= maxRoiEccentricity); %  Enforce Elliptical SHAPE
% 	bwRP = bwRP([bwRP.Perimeter]./([bwRP.Area]) < maxPerimOverSqArea); %  Enforce LOOSELY CIRCULAR/SQUARE SHAPE
% 	bwRP = bwRP([bwRP.Perimeter]./sqrt([bwRP.Area]) > minPerimOverSqArea); %  Enforce NON-HOLINESS (SELF-FULFILLMENT?)
    bwRP = bwRP([bwRP.Orientation] < 10 & [bwRP.Orientation] > -65); % orientation of cells
	if isempty(bwRP)
		continue
	end
	% FILL 'REGIONOFINTEREST' CLASS OBJECTS (several per frame)
    % Using method from Class Object
	frameROI{k,1} = RegionOfInterest(bwRP);
	set(frameROI{k,1},...
		'FrameIdx',frameNum(k),...
		'FrameSize',frameSize);
end
roi = cat(1,frameROI{:});
% if nargin<2
% 	info = [];
% 	frameNum = 1:N;
% else
% 	frameNum = cat(1,info.frame);
% end
% frameROI = cell(N,1);
% parfor kp = 1:N
% 	% EVALUATE CONNECTED COMPONENTS FROM BINARY 'BLOBBED' MATRIX: ENFORCE MORPHOLOGY RESTRICTIONS
% 	bwRP =  regionprops(bwmask(:,:,kp),...
% 		'Centroid', 'Area', 'BoundingBox','PixelIdxList',...
% 		'Eccentricity', 'Orientation','Perimeter');
% 	bwRP = bwRP([bwRP.Area] >= minRoiPixArea);	%	Enforce MINIMUM SIZE
%  	bwRP = bwRP([bwRP.Area] <= maxRoiPixArea);	%	Enforce MAXIMUM SIZE
% 	bwRP = bwRP([bwRP.Eccentricity] >= maxRoiEccentricity); %  Enforce PLUMP SHAPE
% % 	bwRP = bwRP([bwRP.Perimeter]./sqrt([bwRP.Area]) < maxPerimOverSqArea); %  Enforce LOOSELY CIRCULAR/SQUARE SHAPE
% 	bwRP = bwRP([bwRP.Perimeter]./sqrt([bwRP.Area]) > minPerimOverSqArea); %  Enforce NON-HOLINESS (SELF-FULFILLMENT?)
% %     bwRP = bwRP([bwRP.Orientation] > 0 & [bwRP.Orientation] < 80);
% 	if isempty(bwRP)
% 		continue
% 	end
% 	% FILL 'REGIONOFINTEREST' CLASS OBJECTS (several per frame)
% 	frameROI{kp,1} = RegionOfInterest(bwRP);
% 	set(frameROI{kp,1},...
% 		'FrameIdx',frameNum(kp),...
% 		'FrameSize',frameSize);
% end
% roi = cat(1,frameROI{:});
% ------------ SUBFUNCTIONS -------------------
% FUNCTION TO MAKE BINARY MASK WITH ADAPTIVE THRESHOLD
% 	function [bw, sigThresh]  = getAdaptiveHotspots(diffImage, sigThresh)
% 		coverageMaxRatio = .025; %  .01 = 10K pixels (15-30 cells?)
% 		coverageMinPixels = 300; % previous values: 500, 250
% 		thresholdStep = 1;
% 		% PREVENT RACE CONDITION OR NON-SETTLING THRESHOLD
% 		persistent depth
% 		if isempty(depth)
% 			depth = 0;
% 		else
% 			depth = depth + 1;
% 			% 		 thresholdStep = 1 + depth;
% 		end
% 		recursionLim = 250;
% 		if depth > recursionLim
% 			warning('Recursion limit exceeded')
% 			depth = 0;
% 			bw = false(size(diffImage));
% 			sigThresh = gpuArray( stat.Min + uint16( stat.Std.*stdOverMin ));% NEW, (reset)
% 			return
% 		end
% 		% USE THRESHOLD MATRIX TO MAKE BINARY IMAGE, THEN APPLY MORPHOLOGICAL OPERATIONS
% 		diffImage = gpuArray(diffImage);
% 		bw = diffImage > sigThresh;
% 		% changed from: bw = imclose(imopen( bw, S.disk6), S.disk4);
% 		bw = gather(bwmorph(bwmorph(bwmorph( bw, 'open'), 'shrink'), 'majority'));%4
% 		% can also try: 'hbreak'  'shrink' 'fill'  'open' gpuArray
% 		% CHECK FOR OVER/UNDER-THRESHOLDING
% % 		numPix = numel(bw);
% % 		sigThreshPix = sum(bw(:));
% % 		binaryCoverage = sigThreshPix/numPix;
% % 		if binaryCoverage > coverageMaxRatio
% % 			sigThresh = sigThresh + thresholdStep;
% % 			[bw, sigThresh]  = getAdaptiveHotspots(diffImage, sigThresh);
% % 		elseif sigThreshPix < coverageMinPixels
% % 			sigThresh = sigThresh - thresholdStep;
% % 			[bw, sigThresh]  = getAdaptiveHotspots(diffImage, sigThresh);
% % 		else
% % 			depth = 0;
% % 		end
% 	end

end