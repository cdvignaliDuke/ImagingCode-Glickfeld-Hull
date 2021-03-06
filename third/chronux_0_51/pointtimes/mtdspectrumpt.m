function [dS,f]=mtdspectrumpt(data,phi,params,t)
% Multi-taper spectral derivative - point process times
%
% Usage:
%
% [dS,f]=mtdspectrumpt(data,phi,params,t)
% Input: 
%   Note that all times can be in arbitrary units. But the units have to be
%   consistent. So, if E is in secs, win, t have to be in secs, and Fs has to
%   be Hz. If E is in samples, so are win and t, and Fs=1. In case of spike
%   times, the units have to be consistent with the units of data as well.
%       data        (structure array of spike times with dimension channels/trials; also accepts 1d array of spike times) -- required
%       phi         (angle for evaluation of derivative) -- required.
%                       e.g. phi=[0,pi/2] giving the time and frequency derivatives
%       params: structure with fields tapers, pad, Fs, fpass, trialave
%       -optional
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
%           trialave (average over trials when 1, don't average when 0) -
%           optional. Default 0
%       t        (time grid over which the tapers are to be calculated:
%                      this argument is useful when calling the spectrum
%                      calculation routine from a moving window spectrogram
%                      calculation routine). If left empty, the spike times
%                      are used to define the grid.
% Output:
%       dS      (spectral derivative in form phi x frequency x channels/trials if trialave=0; function of phi x frequency if trialave=1)
%       f       (frequencies)
if nargin < 2; error('Need data and angle'); end;
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
dt=1/Fs; % sampling time
if nargin < 4;
   [mintime,maxtime]=minmaxsptimes(data);
   t=mintime:dt:maxtime+dt; % time grid for prolates
end;
N=length(t); % number of points in grid for dpss
nfft=2^(nextpow2(N)+pad); % number of points in fft of prolates
[f,findx]=getfgrid(Fs,nfft,fpass); % get frequency grid for evaluation
tapers=dpsschk(tapers,N,Fs); % check tapers
K=size(tapers,2);
[J,Msp,Nsp]=mtfftpt(data,tapers,nfft,Fs,t,f,findx); % mt fft for point process times
S=squeeze(mean(J(:,1:K-1,:).*conj(J(:,2:K,:)),2));
if trialave; S=squeeze(mean(S,2));end;
nphi=length(phi);
for p=1:nphi;
    dS(p,:,:)=real(exp(i*phi(p))*S);
end;
dS=squeeze(dS);
