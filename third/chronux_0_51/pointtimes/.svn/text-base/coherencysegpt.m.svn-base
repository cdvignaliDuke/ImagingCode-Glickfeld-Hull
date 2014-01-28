function [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegpt(data1,data2,win,params,segave,fscorr)
% Multi-taper coherency computed by segmenting two univariate point processes
% into chunks
%
% Usage:
% [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencysegpt(data1,data2,win,params,segave,fscorr)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data1 (1d structure array of spike times; also accepts 1d array of spike times) -- required
%       data2 (1d structure array of spike times; also accepts 1d array of spike times) -- required
%       win   (length of segments) - required
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
%       segave - optional 0 for don't average over segments, 1 for average - default
%       fscorr   (finite size corrections, 0 (don't use finite size corrections) or 1 (use finite size corrections) - optional
%                (available only for spikes). Defaults 0.
% Output:
%       C (magnitude of coherency - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       phi (phase of coherency - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S12 (cross spectrum -  frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S1 (spectrum 1 - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S2 (spectrum 2 - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       f (frequencies)
%       zerosp (1 for segments where no spikes were found, 0 otherwise)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phierr (error bars for phi) - only for err(1)>=1
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)


if nargin < 3; error('Need data1 and data2 and size of segment'); end;
if nargin < 4; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 5 || isempty(segave); segave=1;end;
if nargin < 6 || isempty(fscorr); fscorr=0; end;

if nargout > 9 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
if nargout > 7 && err(1)==0;
    error('Errors computed only if err(1) is not equal to zero');
end;

[N1,C1,N2,C2]=check_consistency(data1,data2);
[mintime1,maxtime1]=minmaxsptimes(data1);
[mintime2,maxtime2]=minmaxsptimes(data2);
mintime=min(mintime1,mintime2);
maxtime=max(maxtime1,maxtime2);
dt=1/Fs;
t=mintime:dt:maxtime+dt; % time grid for prolates
N=length(t); % number of points in grid for dpss

E=mintime:win:maxtime; % fictitious event triggers
win=[0 win]; % use window length to define left and right limits of windows around triggers
data1=createdatamatpt(data1,E,win); % segmented data 1
data2=createdatamatpt(data2,E,win); % segmented data 2
params.trialave=segave;
if err(1)==0;
   [C,phi,S12,S1,S2,f,zerosp]=coherencypt(data1,data2,params,fscorr); % compute coherency for segmented data
elseif err(1)==1;
   [C,phi,S12,S1,S2,f,zerosp,confC,phierr]=coherencypt(data1,data2,params,fscorr); % compute coherency for segmented data
elseif err(1)==2;
   [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencypt(data1,data2,params,fscorr); % compute coherency for segmented data
end;
