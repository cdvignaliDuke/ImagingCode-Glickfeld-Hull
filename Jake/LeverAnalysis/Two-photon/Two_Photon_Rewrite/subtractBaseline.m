function [data, lastFrame] = subtractBaseline(data, lastFrame)
if nargin < 2
	lastFrame = [];
end						% Data may be approximately between (.5, 1.75)
fprintf('Subtracting Baseline \n')
% SUBTRACT RESULTING BASELINE THAT STILL EXISTS IN NEUROPIL
% dataRange = range(data(:));
% imFloor =  getNearMin(data) - lowBuf;
% fprintf('\t->Adding %3.3g to input\n', -imFloor)
% data = data - imFloor;
activityImage = imfilter(range(data,3), fspecial('average',201), 'replicate');
npMask = double(activityImage) < median(activityImage(:));
npPixNum = sum(npMask(:));
npBaseline = sum(sum(bsxfun(@times, data, cast(npMask,'like',data)), 1), 2) ./ npPixNum; %average of pixels in mask
cellMask = ~npMask;
cellPixNum = sum(cellMask(:));
cellBaseline = sum(sum(bsxfun(@times, data, cast(cellMask,'like',data)), 1), 2) ./ cellPixNum;
% npBaseline = npBaseline(:);
baselineOffset = log(mean(cellBaseline(:))+1);
data = cast( exp( bsxfun(@minus,...
	log(single(data)+1) + baselineOffset,...
	log(single(npBaseline)+1))) - 1, 'like', data) - 1;
% REMOVE BASELINE SHIFTS BETWEEN FRAMES (TODO: untested, maybe move to subtractBaseline)
if ~isempty(lastFrame)
	firstFrame = median(data(:,:,1:8), 3);
	interFileDif = single(firstFrame) - single(lastFrame);
	fileRange = range(data,3);
	baselineShift = double(mode(interFileDif(fileRange < median(fileRange(:)))));
	fprintf('\t->Applying baseline-shift: %3.3g\n',-baselineShift)
	data = data - baselineShift;
end
lastFrame = median(data(:,:,end-7:end), 3);

end