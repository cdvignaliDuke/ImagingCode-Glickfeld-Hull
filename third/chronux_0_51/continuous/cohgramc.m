function [C,phi,S12,S1,S2,t,f,confC,phierr,Cerr]=cohgramc(data1,data2,movingwin,params)
% Multi-taper time-frequency coherence,cross-spectrum and individual spectra - continuous processes
%
% Usage:
%
% [C,phi,S12,S1,S2,t,f,confC,phierr,Cerr]=cohgramc(data1,data2,movingwin,params)
% Input: 
% Note units have to be consistent. Thus, if movingwin is in seconds, Fs
% has to be in Hz. see chronux.m for more information.
%
%       data1 (in form samples x trials) -- required
%       data2 (in form samples x trials) -- required
%       movingwin (in the form [window winstep] -- required
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
% Output:
%       C (magnitude of coherency time x frequencies x trials for trialave=0; time x frequency for trialave=1)
%       phi (phase of coherency time x frequencies x trials for no trial averaging; time x frequency for trialave=1)
%       S12 (cross spectrum - time x frequencies x trials for no trial averaging; time x frequency for trialave=1)
%       S1 (spectrum 1 - time x frequencies x trials for no trial averaging; time x frequency for trialave=1)
%       S2 (spectrum 2 - time x frequencies x trials for no trial averaging; time x frequency for trialave=1)
%       t (time)
%       f (frequencies)
%       confC (confidence level for C at 1-p %) - only for err(1)>=1
%       phierr (error bars for phi) - only for err(1)>=1
%       Cerr  (Jackknife error bars for C - use only for Jackknife - err(1)=2)

if nargin < 3; error('Need data1 and data2 and window parameters'); end;
if nargin < 4; params=[];end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargout > 9 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
if nargout > 7 && err(1)==0;
%   Errors computed only if err(1) is nonzero. Need to change params and run again.
    error('When errors are desired, err(1) has to be non-zero.');
end;
[N1,C1,N2,C2]=check_consistency(data1,data2);
N=N1;

Nwin=round(Fs*movingwin(1)); % number of samples in window
Nstep=round(movingwin(2)*Fs); % number of samples to step through
nfft=2^(nextpow2(Nwin)+pad);
f=getfgrid(Fs,nfft,fpass); 
params.tapers=dpsschk(tapers,Nwin,Fs); % check tapers

winstart=1:Nstep:N-Nwin+1;
nw=length(winstart);
for n=1:nw;
   indx=winstart(n):winstart(n)+Nwin-1;
   datawin1=data1(indx,:);datawin2=data2(indx,:);
   if nargout==10;
     [c,ph,s12,s1,s2,f,confc,phie,cerr]=coherencyc(datawin1,datawin2,params);
     confC=confc;
     phierr(1,n,:,:)=squeeze(phie(1,:,:));
     phierr(2,n,:,:)=squeeze(phie(2,:,:));
     Cerr(1,n,:,:)=squeeze(cerr(1,:,:));
     Cerr(2,n,:,:)=squeeze(cerr(2,:,:));
   elseif nargout==9;
     [c,ph,s12,s1,s2,f,confc,phie]=coherencyc(datawin1,datawin2,params);
     confC=confc;
     phierr(1,n,:,:)=squeeze(phie(1,:,:));
     phierr(2,n,:,:)=squeeze(phie(2,:,:));
   else
     [c,ph,s12,s1,s2,f]=coherencyc(datawin1,datawin2,params);
   end;
   C(n,:,:)=c;
   S12(n,:,:)=s12;
   S1(n,:,:)=s1;
   S2(n,:,:)=s2;
   phi(n,:,:)=ph;
end;
C=squeeze(C); phi=squeeze(phi);S12=squeeze(S12); S1=squeeze(S1); S2=squeeze(S2);
if nargout==10;Cerr=squeeze(Cerr);end;
if nargout>=9; phierr=squeeze(phierr);end
winmid=winstart+round(Nwin/2);
t=winmid/Fs;
