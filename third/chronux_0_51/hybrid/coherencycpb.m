function [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpb(data1,data2,params,fscorr)
% Multi-taper coherency,cross-spectrum and individual spectra - continuous and binned point process data
%
% Usage:
%
% [C,phi,S12,S1,S2,f,zerosp,confC,phierr,Cerr]=coherencycpb(data1,data2,params,fscorr)
% Inputs: 
%       data1 (continuous data in form samples x trials) -- required
%       data2 (binned spike data in form samples x trials) -- required
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
% Outputs:
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

if nargin < 2; error('Need data1 and data2'); end
if nargin < 3; params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin < 4 || isempty(fscorr); fscorr=0; end;

if nargout > 7 && err(1)==0;
%Cannot compute error bars if err(1)=0. Need to change params and run again. 
    error('When errors are desired, err(1) has to be non-zero.');    
end;
if nargout > 9 && err(1)~=2; 
    error('Cerr computed only for Jackknife. Correct inputs and run again');
end;
[N1,C1,N2,C2]=check_consistency(data1,data2);
N=N1;
zerosp=zeros(1,C1); % intialize the zerosp variable

nfft=2^(nextpow2(N)+pad);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers
J1=mtfftc(data1,tapers,nfft,Fs); % fourier transform of continuous data
[J2,Msp2,Nsp2]=mtfftpb(data2,tapers,nfft); % fourier transform of point process data
zerosp(find(Nsp2==0))=1; % set the trials where no spikes were found to have zerosp=1
J1=J1(findx,:,:);
J2=J2(findx,:,:);
S12=squeeze(mean(conj(J1).*J2,2));
S1=squeeze(mean(conj(J1).*J1,2));
S2=squeeze(mean(conj(J2).*J2,2));
if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;
C12=S12./sqrt(S1.*S2);
C=abs(C12);
phi=angle(C12);
if nargout==10; 
  if fscorr==1; 
     [confC,phierr,Cerr]=coherr(C,J1,J2,err,trialave,[],Nsp2);
  else
     [confC,phierr,Cerr]=coherr(C,J1,J2,err,trialave);
  end;
elseif nargout==9;
  if fscorr==1; 
     [confC,phierr]=coherr(C,J1,J2,err,trialave,[],Nsp2);
  else
     [confC,phierr]=coherr(C,J1,J2,err,trialave);
  end;
end;
