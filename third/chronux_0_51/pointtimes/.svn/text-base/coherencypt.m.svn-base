function [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencypt(data1,data2,params,fscorr,t)
% Multi-taper coherency - point process times
%
% Usage:
%
% [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencypt(data1,data2,params,fscorr,t)
% Input: 
%       data1  (structure array of spike times with dimension trials; also accepts 1d array of spike times) -- required
%       data2  (structure array of spike times with dimension trials; also accepts 1d array of spike times) -- required
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
%           trialave (average over trials when 1, don't average when 0) - optional. Default 0
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
%       t        (time grid over which the tapers are to be calculated:
%                      this argument is useful when calling the spectrum
%                      calculation routine from a moving window spectrogram
%                      calculation routine). If left empty, the spike times
%                      are used to define the grid.
% Output:
%       C (magnitude of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       phi (phase of coherency - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S12 (cross spectrum -  frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S1 (spectrum 1 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       S2 (spectrum 2 - frequencies x trials if trialave=0; dimension frequencies if trialave=1)
%       f (frequencies)
%       zerosp (1 for trials where no spikes were found, 0 otherwise)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phierr (error bars for phi) - only for err(1)>=1
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)
if nargin < 2; error('Need data1 and data2'); end;
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 4 || isempty(fscorr); fscorr=0; end;
if nargin < 5 || isempty(t); 
   [mintime1,maxtime1]=minmaxsptimes(data1);
   [mintime2,maxtime2]=minmaxsptimes(data2);
   mintime=min(mintime1,mintime2);
   maxtime=max(maxtime1,maxtime2);
   dt=1/Fs;
   t=mintime:dt:maxtime+dt; % time grid for prolates
end;
   
if nargout > 9 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs or outputs and run again');
end;
if nargout > 7 && err(1)==0;
    error('Errors computed only if err(1) is not equal to zero');
end;

[N1,C1,N2,C2]=check_consistency(data1,data2);

N=length(t); % number of points in grid for dpss
nfft=2^(nextpow2(N)+pad); % number of points in fft of prolates
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
[J1,Msp1,Nsp1]=mtfftpt(data1,tapers,nfft,Fs,t,f,findx);
[J2,Msp2,Nsp2]=mtfftpt(data2,tapers,nfft,Fs,t,f,findx);
zerosp=zeros(1,C1); % initialize the zerosp variable
zerosp(Nsp1==0 | Nsp2==0)=1; % set the zerosp variable
S12=squeeze(mean(conj(J1).*J2,2));
S1=squeeze(mean(conj(J1).*J1,2));
S2=squeeze(mean(conj(J2).*J2,2));
if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;
C12=S12./sqrt(S1.*S2);
C=abs(C12);
phi=angle(C12);
if nargout==10; 
    if fscorr==1; 
       [confC,phierr,Cerr]=coherr(C,J1,J2,err,trialave,Nsp1,Nsp2);
    else
       [confC,phierr,Cerr]=coherr(C,J1,J2,err,trialave);
    end;
elseif nargout==9;
    if fscorr==1; 
        [confC,phierr]=coherr(C,J1,J2,err,trialave,Nsp1,Nsp2);
    else
        [confC,phierr]=coherr(C,J1,J2,err,trialave);
    end;
end;
