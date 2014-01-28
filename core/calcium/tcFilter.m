function tc = tcFilter(tcIn, varargin)
%TCFILTER (calcium): on tc, do high-pass, dF/F, trial averaging, final smooth
%   TC = TCFILTER(TCIN, STROPTS)
%
%   tc is a var of size [nTimepts,nCells]
%
%   Run tcFilter without arguments to list the optional and default arguments
% 
%   See also SMOOTH, FILTER, BUTTER, INTERP
%
%$Id: tcFilter.m 159 2008-04-01 20:39:32Z vincent $

uDefs = { 'FrameTimeMs', [], ...
          'DoHighpassFilter', true, ...
          'HighpassHz', 0.05, ...
          'DoHighpassRestoreMean', true, ...
          'FilterType', 'firls', ...
          'FirlsFilterOrder', 300, ...
          'BaselineFrames', [], ...
          'FilterInterpFrames', [], ...
          'DoDFOF', true, ...
          'FrameNsToAverage', [], ...
          'DoPostSmooth', false, ...   % see SMOOTH.m
          'SmoothMethod', 'lowess', ...
          'SmoothSpan', 5, ...
        };

uo = stropt2struct(stropt_defaults(uDefs, varargin));

if nargin ==0
    display(uo)
    return;
end


%%

[nTimepts,nCells] = size(tcIn);
if nTimepts == 1, error('tcIn must be size [nTimepts,nCells]'); end

tc = double(tcIn);

% precompute some things
if uo.DoDFOF || (uo.DoHighpassFilter && uo.DoHighpassRestoreMean)
    if isempty(uo.BaselineFrames)
        error('Must specify baseline frames');
    end
    baseVals = mean(tc(uo.BaselineFrames,:),1);
    baseCols = repmat(baseVals, [nTimepts,1]);
end

%% filtering
if uo.DoHighpassFilter
    if isempty(uo.FrameTimeMs)
        error('Must specifiy FrameTimeMs in order to filter');
    end
    if ~strcmp(uo.FilterType, 'interp') && ~isempty(uo.FilterInterpFrames)
        error('Cannot specify FilterInterpFrames for non-''interp'' filters');
    end


    switch uo.FilterType
      case 'butter'
        % butterworth
        Fs = 1000/uo.FrameTimeMs;
        Fstop = uo.HighpassHz / (Fs/2);
        N = 6;
        [b,a] = butter(N, [Fstop], 'high');
        tcF = filtfilt(b,a,tc);
      case 'firls'
        % finite impulse reponse
        Fs = 1000/uo.FrameTimeMs;
        fpts = [0 uo.HighpassHz/2 uo.HighpassHz Fs/2] / (Fs/2);
        apts = [0 0 1 1];
        wpts = [1 1];
        N = uo.FirlsFilterOrder;
        [b,a] = firls(N, fpts, apts, wpts);
        tcF = filtfilt(b,a,tc);
      case 'interp'
        % compute filtered baseline
        allLowF = interpLowpass(tc, uo.FilterInterpFrames, ...
                                'butter', 1000/uo.FrameTimeMs, ...
                                uo.HighpassHz);
        % subtract this baseline
        tcF = tc - allLowF;
      case 'fft' % fft-based method by Kenichi Ohki
          % cutoff period in frames 
          cutoffperiod = floor(1/uo.HighpassHz/(uo.FrameTimeMs*1e-3))
          tcF = tcLowCut (tc, cutoffperiod,'gaussian')
            
      otherwise 
        error('unknown filter type %s', uo.FilterType);
    end

    if uo.DoHighpassRestoreMean
        % add mean back in
        hpMean = mean(tcF(uo.BaselineFrames,:),1);
        tcF = tcF + baseCols - repmat(hpMean, [nTimepts,1]);
    end

    tc = tcF;
end

%% dF/F
if uo.DoDFOF
    tc = (tc - baseCols) ./ baseCols;
end

%% averaging
if ~isempty(uo.FrameNsToAverage)
    [nFramesInAverage nOutFr] = size(uo.FrameNsToAverage);
    avgTc = repmat(NaN, [nOutFr,nCells]);
    for iF = 1:nOutFr
        avgTc(iF,:) = mean(tc(uo.FrameNsToAverage(:,iF),:),1);
    end
    tc = avgTc;
end

%% post-average smoothing
if uo.DoPostSmooth
    for iC = 1:nCells
        tc(:,iC) = smooth(tc(:,iC), uo.SmoothSpan, uo.SmoothMethod);
    end
end

% old effort to interp filter: subset, butterworth, re-interpolate.
% didn't work because interpolation is hard
%
% $$$       case 'interp'
% $$$ 
% $$$         % select a subset of frames, lowpass, interpolate, subtract
% $$$         assert(~isempty(uo.FilterInterpFrames));
% $$$         assert(all(diff(uo.FilterInterpFrames) > 0), ...
% $$$                'FilterInterpFrames must be monotonically increasing');
% $$$         
% $$$ 
% $$$         tcPre = tc(uo.FilterInterpFrames);
% $$$         
% $$$         % lowpass, butter
% $$$         Fs = 1000/uo.FrameTimeMs;
% $$$         Fstop = uo.HighpassHz / (Fs/2);
% $$$         N = 6;
% $$$         [b,a] = butter(N, [Fstop], 'low');
% $$$         lowF = filtfilt(b,a,tcPre')';
% $$$ 
% $$$         % interp to all frames
% $$$         allLowF = interp1(uo.FilterInterpFrames, lowF, 1:nTimepts, 'pchip');
% $$$         keyboard
% $$$         % subtract
% $$$         tcF = tc - allLowF;
% $$$         
% $$$         %% test, interp first
% $$$         r = interp1(uo.FilterInterpFrames, tc(:,uo.FilterInterpFrames), ...
% $$$                     1:nTimepts, 'pchip');
% $$$         
% $$$         %% end
