function [S,t,f,R,Serr]=mtspecgrampt(data,movingwin,params,fscorr)
% Multi-taper time-frequency spectrum - point process times
%
% Usage:
%
% [S,t,f,R,Serr]=mtspecgrampt(data,movingwin,params,fscorr)
% Input: 
%       data        (structure array of spike times with dimension channels/trials; also accepts 1d array of spike times) -- required
%       movingwin         (in the form [window,winstep] i.e length of moving
%                                                 window and step size.
%                                                 
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
%           trialave (average over trials/channels when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
%
% Output:
%       S       (spectrogram with dimensions time x frequency x channels/trials if trialave=0; dimensions time x frequency if trialave=1)
%       t       (times)
%       f       (frequencies)
%
%       Serr    (error bars) - only if err(1)>=1

if nargin < 2; error('Need data and window parameters'); end;
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 4 || isempty(fscorr); fscorr=0; end;
if nargout > 4 && err(1)==0; error('Cannot compute errors with err(1)=0'); end;

[mintime,maxtime]=minmaxsptimes(data);
tn=(mintime+movingwin(1)/2:movingwin(2):maxtime-movingwin(1)/2);
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=2^(nextpow2(Nwin)+pad);
[f,findx]=getfgrid(Fs,nfft,fpass); 
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
nw=length(tn);
for n=1:nw;
   t=linspace(tn(n)-movingwin(1)/2,tn(n)+movingwin(1)/2,Nwin);
   datawin=extractdatapt(data,[t(1) t(end)]);
   if nargout==5;
     [s,f,r,serr]=mtspectrumpt(datawin,params,fscorr,t);
     Serr(1,n,:,:)=squeeze(serr(1,:,:));
     Serr(2,n,:,:)=squeeze(serr(2,:,:));
   else
     [s,f,r]=mtspectrumpt(datawin,params,fscorr,t);
   end;
   S(n,:,:)=s;
   R(n,:)=r;
end;
t=tn;
S=squeeze(S); R=squeeze(R); if nargout==5; Serr=squeeze(Serr);end
