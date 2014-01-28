function [S,t,f,R,Serr]=mtspecgrampb(data,movingwin,params,fscorr)
% Multi-taper time-frequency spectrum - binned point process
%
% Usage:
%
% [S,t,f,R,Serr]=mtspecgrampb(data,movingwin,params,fscorr)
% Input: 
%       data        (in form samples x channels/trials) -- required
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
%           trialave (average over trials/channnels when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       S       (spectrum in form time x frequency x channels/trials for trialave=0; or as a function of frequency if trialave=1)
%       t       (times)
%       f       (frequencies)
%       R       (rate)
%       Serr    (error bars) - only for err(1)>=1

if nargin < 2; error('Need data and window parameters'); end;
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 4 || isempty(fscorr); fscorr=0; end;
if nargout > 4 && err(1)==0; 
%    error('Cannot compute errors with err(1)=0'); 
     error('When Serr is desired, err(1) has to be non-zero.');
end;
N=size(data,1);
Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=2^(nextpow2(Nwin)+pad);
[f,findx]=getfgrid(Fs,nfft,fpass); 
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers
whos data
winstart=1:Nstep:N-Nwin+1;
nw=length(winstart);
for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin=data(indx,:);
   if nargout==5;
     [s,f,r,serr]=mtspectrumpb(datawin,params,fscorr);
     Serr(1,n,:,:)=squeeze(serr(1,:,:));
     Serr(2,n,:,:)=squeeze(serr(2,:,:));
   else
     [s,f,r]=mtspectrumpb(datawin,params,fscorr);
   end;
   S(n,:,:)=s;
   R(n,:)=r';
end;
winmid=winstart+round(Nwin/2);
t=winmid/Fs;
S=squeeze(S); R=squeeze(R); if nargout==5; Serr=squeeze(Serr);end
