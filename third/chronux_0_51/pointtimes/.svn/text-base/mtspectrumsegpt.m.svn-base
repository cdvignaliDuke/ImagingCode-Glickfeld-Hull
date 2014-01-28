function [S,f,R,varS,zerosp,C,Serr]=mtspectrumsegpt(data,win,params,segave,fscorr)
% Multi-taper segmented spectrum for a univariate binned point process
%
% Usage:
%
% [S,f,R,varS,zerosp,C,Serr]=mtspectrumsegpt(data,win,params,segave,fscorr)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data (structure array of one channel of spike times; also accepts 1d column vector of spike times) -- required
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
%       segave - (0 for don't average over segments, 1 for average) - optional - default  1
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       S       (spectrum in form frequency x segments if segave=0; function of frequency if segave=1)
%       f       (frequencies)
%       R       (spike rate)
%       varS    (variance of the spectrum as a function of frequency)
%       zerosp  (0 for segments in which spikes were found, 1 for segments
%       C       (covariance matrix of the log spectrum - frequency x
%       frequency matrix)
%       Serr    (error bars) - only if err(1)>=1

if nargin < 2; error('Need data and segment information'); end;
if nargin < 3; params=[]; end;
if nargin < 4 || isempty(segave); segave=1; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 5 || isempty(fscorr); fscorr=0;end;

if nargout > 4 && err(1)==0; error('cannot compute error bars with err(1)=0; change params and run again'); end;

if isstruct(data);
   fnames=fieldnames(data);
   eval(['dtmp=data.' fnames{1} ';'])
else
   dtmp=data(:);
end;
T=max(dtmp); % total length of data - USED????
minT=min(dtmp); 
dt=1/Fs; % sampling interval
E=minT:win:T-win; % fictitious event triggers
win=[0 win]; % use window length to define left and right limits of windows around triggers
dtmp=createdatamatpt(dtmp,E,win); % create segmented data set
zerosp=zeros(1,size(data,1));
[mintime,maxtime]=minmaxsptimes(dtmp);
dt=1/Fs; % sampling time
t=mintime-dt:dt:maxtime+dt; % time grid for prolates
N=length(t); % number of points in grid for dpss
nfft=2^(nextpow2(N)+pad); % number of points in fft of prolates
[f,findx]=getfgrid(Fs,nfft,fpass); % get frequency grid for evaluation
tapers=dpsschk(tapers,N,Fs); % check tapers
[J,Msp,Nsp]=mtfftpt(dtmp,tapers,nfft,Fs,t,f,findx);% mt fft for point process times
zerosp(find(Nsp==0))=1;
R=Msp*Fs;
S=squeeze(mean(conj(J).*J,2)); % spectra of non-overlapping segments (averaged over tapers)
lS=log(S); % log spectrum
C=cov(lS')'; % covariance matrix of the log spectrum
varS=diag(C); % variance of fthe log spectrum
if segave==1; S=squeeze(mean(S,2)); R=squeeze(mean(R)); end; % mean of the spectrum averaged across segments and tapers
if nargout==6; 
   if fscorr==1;
      Serr=specerr(S,J,err,trialave,Nsp);
   else
      Serr=specerr(S,J,err,trialave);
   end;
end;
