function [S,f,R,Serr]=mtspectrumtrigpb(data,E,win,params,fscorr)
% Multi-taper event triggered time-frequency spectrum - binned point process
%
% Usage:
%
% [S,f,R,Serr]=mtspectrumtrigpb(data,E,win,params,fscorr)
% Input: 
%       data        (single channel data in a column vector) -- required
%       E           (event times) - required
%       win         (in the form [winl winr] i.e window around each event
%                                                 required
%                                                 Note that units here have
%                                                 to be consistent with
%                                                 units of Fs
%       params: structure with fields tapers, pad, Fs, fpass, err, trialave
%       - optional
%           tapers (precalculated tapers from dpss, or in the form [NW K] e.g [3 5]) -- optional. If not 
%                                                 specified, use [NW K]=[3 5]
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	 e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	 to 512 points; if PAD = 2, we pad the FFT
%			      	 to 2048 points, etc.
%           Fs   (sampling frequency) - optional. Default 1.
%           fpass    (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%           err  (error calculation [1 p] - Theoretical error bars; [2 p] - Jackknife error bars
%                                   [0 p] or 0 - no error bars) - optional. Default 0.
%           trialave (average over events when 1, don't average when 0) -
%           optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       S       (triggered spectrum in form frequency x events for trialave=0, or as a function of frequency for trialave=1)
%       f       (frequencies)
%       R       (spike rate)
%       Serr    (error bars) - only for err(1)>=1

if nargin < 3; error('Need data, events and window parameters'); end;
if nargin < 2; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 5 || isempty(fscorr); fscorr=0; end;
if nargout > 3 && err(1)==0; 
%   Cannot compute errors if err(1)=0. Need to change params and run again. 
    error('When Serr is desired, err(1) has to be non-zero.');
end;

data=createdatamatpb(data,E,Fs,win); 
if nargout==4; [S,f,R,Serr]=mtspectrumpb(data,params,fscorr);
else [S,f,R]=mtspectrumpb(data,params,fscorr);end
