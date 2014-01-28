function [C,phi,S12,S1,S2,f,confC,phierr,Cerr]=coherencysegc(data1,data2,win,params)
% Multi-taper coherency,cross-spectrum and individual spectra computed by segmenting two univariate time series into chunks - continuous process
%
% Usage:
% [C,phi,S12,S1,S2,f,confC,phierr,Cerr]=coherencysegc(data1,data2,win,params)
% Input: 
% Note units have to be consistent. See chronux.m for more information.
%       data1 (column vector) -- required
%       data2 (column vector) -- required
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
% Output:
%       C (magnitude of coherency - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       phi (phase of coherency - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S12 (cross spectrum -  frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S1 (spectrum 1 - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       S2 (spectrum 2 - frequencies x segments if segave=0; dimension frequencies if segave=1)
%       f (frequencies)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phierr (error bars for phi) - only for err(1)>=1
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)

if nargin < 3; error('Need data1 and data2 and size of segment'); end;
if nargin < 4; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargout > 8 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
if nargout > 6 && err(1)==0;
%   Errors computed only if err(1) is nonzero. Need to change params and run again. 
    error('When errors are desired, err(1) has to be non-zero.');
end;
if size(data1,2)~=1 || size(data2,2)~=1; error('works for only univariate time series'); end;

[N1,C1,N2,C2]=check_consistency(data1,data2);
N=N1;

dt=1/Fs; % sampling interval
T=N*dt; % length of data in seconds
E=0:win:T-win; % fictitious event triggers
win=[0 win]; % use window length to define left and right limits of windows around triggers
data1=createdatamatc(data1,E,Fs,win); % segmented data 1
data2=createdatamatc(data2,E,Fs,win); % segmented data 2
params.trialave=1;
params.trialave=1;
if err==0;
   [C,phi,S12,S1,S2,f]=coherencyc(data1,data2,params); % compute coherency for segmented data
elseif err(1)==1;
   [C,phi,S12,S1,S2,f,confC,phierr]=coherencyc(data1,data2,params); % compute coherency for segmented data
elseif err(1)==2;
   [C,phi,S12,S1,S2,f,confC,phierr,Cerr]=coherencyc(data1,data2,params); % compute coherency for segmented data
end;