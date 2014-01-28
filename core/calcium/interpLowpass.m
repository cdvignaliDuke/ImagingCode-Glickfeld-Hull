function [allLowF] = interpLowpass(tcIn, interpFrameNs, filterType, ...
                                origSamplingFreq, lowPassHz)
%INTERP_LOWPASS (stacks-mh): filter a timecourse using only a subset of frames
%
%$Id: interpLowpass.m 164 2008-04-07 04:32:20Z histed $

if nargin<3, filterType = 'butter'; end

[nTimepts,nCells] = size(tcIn);
          
%% check inputs
assert(~isempty(interpFrameNs));
assert(all(diff(interpFrameNs) > 0), ...
       'interpFrameNs must be monotonically increasing');

%% subsample input, and pad on either end to prevent edge effects
tcPre = tcIn(interpFrameNs,:);
allNs = 1:nTimepts;
preNs = allNs(interpFrameNs)';

% pad with mean on either side
nIFrs = length(interpFrameNs);
nPad = min(50, nIFrs-1);
assert(nIFrs > nPad, 'too few timepts');
tcPre = cat(1, ...
            repmat(mean(tcPre(1:nPad,:),1), [nPad 1]), ...
            tcPre, ...
            repmat(mean(tcPre(end-nPad:end,:),1), [nPad 1]));
preNs = cat(1, 1+(-nPad:-1)', preNs, nTimepts+(1:nPad)');

%% do lowpass filter
switch filterType
  case 'butter'
    Fs = origSamplingFreq;
    subSampFract = length(interpFrameNs)/nTimepts;
    % decrease Fs to account for us taking only subsample of frames
    % Fs = Fs ./ subSampFract;   
    %  remove - it makes us too sensitive to fast changes
    Fstop = lowPassHz / (Fs/2);
    N = 6;
    [b,a] = butter(N, [Fstop], 'low');
    lowF = filtfilt(b,a,tcPre);
  otherwise
    error('bad filterType: %s', filterType);
end

%% interp back to all frames
allLowF = interp1(preNs, lowF, 1:nTimepts, 'spline');

if isvector(allLowF)
    allLowF = colvect(allLowF);
end

%keyboard


        
          
