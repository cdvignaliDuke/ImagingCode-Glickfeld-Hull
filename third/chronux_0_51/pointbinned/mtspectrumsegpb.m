function [S,f,R,varS,zerosp,C,Serr]=mtspectrumsegpb(data,win,params,segave,fscorr)
% Multi-taper segmented spectrum for a univariate binned point process
%
% Usage:
%
% [S,f,R,varS,zerosp,C,Serr]=mtspectrumsegpb(data,win,params,segave,fscorr)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (in form samples x 1) -- required
%       win  (duration of the segments) - required. 
%       params: structure with fields tapers, pad, Fs, fpass, err
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
%       segave (1 for averaging across segments, 0 otherwise; default 1)
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       S       (spectrum in form frequency x segments if segave=0; as a function of frequency if segave=1)
%       f       (frequencies)
%       R       (spike rate)
%       varS    (variance of the log spectrum)
%       zerosp  (0 for segments in which spikes were found, 1 for segments
%       in which there are no spikes)
%       C       (covariance matrix of the log spectrum - frequency x
%       frequency matrix)
%       Serr    (error bars) - only for err(1)>=1


if nargin < 2; error('Need data and segment information'); end;
if nargin < 3; params=[]; end;
if nargin < 4 || isempty(segave); segave=1; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 3 || isempty(fscorr); fscorr=0;end;

if nargout > 4 && err(1)==0; 
%   Cannot compute error bars with err(1)=0. Need to change params and run again. 
    error('When Serr is desired, err(1) has to be non-zero.');
end;

N=size(data,1); % total length of data
dt=1/Fs; % sampling interval
T=N*dt; % length of data in seconds
E=0:win:T-win; % fictitious event triggers
win=[0 win]; % use window length to define left and right limits of windows around triggers
data=createdatamatpb(data,E,Fs,win);
zerosp=zeros(1,size(data,2));
N=size(data,1); % length of segmented data
nfft=2^(nextpow2(N)+pad);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
[J,Msp,Nsp]=mtfftpb(data,tapers,nfft);  
zerosp(find(Nsp==0))=1;
J=J(findx,:,:);
R=Msp*Fs;
S=squeeze(mean(conj(J).*J,2)); % spectra of non-overlapping segments (averaged over tapers)
lS=log(S); % log spectrum
C=cov(lS')'; % covariance of the log spectrum - frequency x frequency
varS=diag(C); % variance of the log spectrum
if segave;
  S=squeeze(mean(S,2)); % mean of the spectrum averaged across segments
  R=squeeze(mean(R)); % mean of the spike rate averaged across segments
end
if nargout==6; 
   if fscorr==1;
      Serr=specerr(S,J,err,trialave,Nsp);
   else
      Serr=specerr(S,J,err,trialave);
   end;
end;
