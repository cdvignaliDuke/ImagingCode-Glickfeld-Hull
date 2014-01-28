function [outFrames, preRange, cMap] = make_stim_avgmovie(stack, varargin)
% calculate post-stimulus time movie: 
%   as dF/F or by averaging
%$Id: make_stim_avgmovie.m 297 2008-07-15 21:50:27Z histed $

defs = { 'FrameNsToAverage', [], ...
         'BlockStimNs', [], ...
         'StimLengthFrames', 1, ...
         'MovieType', 'dF', ...  % 'dF', 'dFOF'         
         'DFOFBlockBaseNs', [], ...  %
         'DFOFOutputRange', [-0.5 2.0], ...
         'DFDoLocalContrast', true, ...
         'DoPlot', true, ...
         'PlotColormap', gray(256), ...
         'PlotDoColorbar', false, ...
         'PlotFramesPerSec', 'Default', ...
         'DoSmooth', false, ...
         'SmoothWidth', 2, ...
         'SmoothSigma', 3, ...
         'DecimateTimeBy', [], ...
         'PrintStatus', true };

% $$$          'StimFrPeakOffset', 0, ...
% $$$          'BasePreOffset', 0, ...
% $$$          'PostFrameNs', 1, ...
% $$$          'CLim', [-0.1 0.35], ...
% $$$          'DoMapSmooth', true, ...
% $$$          'MapSmoothWidth', 2, ...
% $$$          'MapSmoothSigma', 3, ...


uo = stropt2struct(stropt_defaults(defs, varargin));
         
%%%%
[nRows nCols nStFrames] = size(stack);
[nStims nAvgFrames] = size(uo.FrameNsToAverage);

if uo.PrintStatus
    fprintf(1, 'Averaging pixels for %d movie frames\n', nAvgFrames);
end
outFrames = repmat(double(0), [nRows, nCols, nAvgFrames]);
for iT = 1:nAvgFrames
    tFrNs = uo.FrameNsToAverage(:,iT);
    switch uo.MovieType
      case { 'dF', 'dFOF' }
        % average frames
        outFrames(:,:,iT) = mean(stack(:,:,tFrNs), 3);
      otherwise
        error('invalid MovieType %s', uo.MovieType);
    end
end

%% post processing
switch uo.MovieType
  case 'dF'
    % local contrast if desired
    if uo.DFDoLocalContrast
        outFrames = stack_localcontrastadj(outFrames, 31, 5);
    end
  case 'dFOF'
    % do dF/F computation
    baseFr = mean(outFrames(:,:,uo.DFOFBlockBaseNs),3);
    baseD = baseFr;
    baseD(baseD<40) = 40;
    outFrames = (outFrames - repmat(baseFr,[1 1 nAvgFrames])) ...
                ./ repmat(baseD,[1 1 nAvgFrames]);  % all prev. doubles
    
end

%% smoothing
if uo.DoSmooth
    for iF = 1:nAvgFrames
        outFrames(:,:,iF) = smooth2(outFrames(:,:,iF), 'gauss', ...
                                    uo.SmoothWidth*[1 1], ...
                                    uo.SmoothSigma.*[1 1]);
    end
end

%% decimate
if ~isempty(uo.DecimateTimeBy)
    nDec = uo.DecimateTimeBy;
    nFrames = size(outFrames,3);
    nRes = floor(nFrames/nDec);
    nToUse = nRes*nDec;
    error('fix this');
    outFrames = mean(reshape(permute(outFrames(:,:,1:nToUse), [3 1 2]), ...
                             [nRes nRows nCols nDec]), ...
                     4);

end

%% normalize
if uo.PrintStatus
    fprintf(1, 'Rescaling movie to 8bits\n');
end
switch uo.MovieType
  case 'dF'
    allMin = min(outFrames(:));
    allMax = max(outFrames(:));
    preRange = [allMin allMax];  % future: add opt to force specified new range
  case 'dFOF'
    preRange = uo.DFOFOutputRange;
end
outFrames = imScale(outFrames, preRange, [1 255], 'uint8');

%% annotate with a square on stim frames
% first colormap; rescale to [1 255]
cMap = interp1([0:(size(uo.PlotColormap,1)-1)], ...
       uo.PlotColormap, 1:255, 'cubic');
cMap = cat(1,[1 1 1],cMap);

nSqPix = round(nRows / 35);
outFrames((end-nSqPix+1):end,1:nSqPix,...
          uo.BlockStimNs:uo.BlockStimNs+uo.StimLengthFrames) = 0;



%% plot it 
if uo.DoPlot
    stackAnimator(outFrames, ...
                  'Colormap', cMap, ...
                  'DoColorbar', true, ...
                  'FramesPerSec', uo.PlotFramesPerSec, ...
                  'CLim', []); % let it figure out best clim from data
end
