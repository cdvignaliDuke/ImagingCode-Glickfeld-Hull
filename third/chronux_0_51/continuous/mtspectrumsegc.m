function [S,f,varS,C,Serr]=mtspectrumsegc(data,win,params,segave)
% Multi-taper segmented spectrum for a univariate continuous process
%
% Usage:
%
% [S,f,varS,C,Serr]=mtspectrumsegc(data,win,params,segave)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (column vector) -- required
%       win  (duration of the segments) - required. 
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
%           trialave - not used
%       segave - optional 0 for don't average over segments, 1 for average - default
%       1
% Output:
%       S       (spectrum in form frequency x segments if segave=0; in the form frequency if segave=1)
%       f       (frequencies)
%       varS    (variance of the log spectrum)
%       C       (covariance matrix of the log spectrum - frequency x
%       frequency matrix)
%       Serr    (error bars) only for err(1)>=1

if nargin < 2; error('Need data and segment information'); end;
if size(data,2)~=1; error('works for only univariate time series'); end;
if nargin < 3 ; params=[]; end;
if nargin < 4 || isempty(segave); segave=1; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargout==4 && err(1)==0; 
%   Errors can't be computed if err(1)=0. Need to change params and run again.
    error('When Serr is desired, err(1) has to be non-zero.');
end;

N=size(data,1); % length of segmented data
dt=1/Fs; % sampling interval
T=N*dt; % length of data in seconds
E=0:win:T-win; % fictitious event triggers
win=[0 win]; % use window length to define left and right limits of windows around triggers
data=createdatamatc(data,E,Fs,win); % segmented data
N=size(data,1); % length of segmented data
nfft=2^(nextpow2(N)+pad);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
J=mtfftc(data,tapers,nfft,Fs); % compute tapered fourier transforms
J=J(findx,:,:); % restrict to specified frequencies
S=squeeze(mean(conj(J).*J,2)); % spectra of non-overlapping segments (average over tapers)
lS=log(S); % log spectrum for non-overlapping segments
C=cov(lS'); % covariance matrix of the log spectrum
varS=diag(C); % variance of the log spectrum as a function of frequency
if segave==1; S=squeeze(mean(S,2)); end; % mean of the spectrum averaged across segments
if nargout==5; 
   Serr=specerr(S,J,err,1);
end;
