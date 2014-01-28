function [sv,sp,fm] = spsvd(data,params,mdkp)
% Space frequency SVD of input data - continuous processes
% Usage: [sv,sp,fm] = spsvd(data,params,mdkp)
% Inputs:
% data       (data matrix in timexchannels form)-required
%       params      structure containing parameters - params has the
%       following fields: tapers, Fs, fpass, pad
%	        tapers 	    (parameters for calculating tapers [NW,K]) - optional. Defaults to [3 5]
%	        Fs 	        (sampling frequency) -- optional. Defaults to 1.
%           fpass       (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	    e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	    to 512 points; if PAD = 2, we pad the FFT
%			      	    to 2048 points, etc.
% mdkp       (number of dimensions to be kept)-optional. Default is the
%               maximum possible modes determined by taper parameters
%
% Outputs:
% sv sp fm  : singular values, space modes, frequency modes


if nargin < 1; error('Need data'); end;
if nargin < 2 || isempty(params); params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);

N=size(data,1);
tapers=dpsschk(tapers,N,Fs);
nfft=2^(nextpow2(N)+pad);% number of points in fft
[N,K]=size(tapers);
if isempty(mdkp) || nargin<5; mdkp=min(K,NCHAN);
elseif mdkp > min(K,NCHAN); error('mdkp has to be less than both K and NCHAN');end;

tvec=1:N;
tvec=tvec*2*pi*i;
f=getfgrid(Fs,nfft,fpass);
nf=length(f);

sp=zeros(nwin,NCHAN,nf,mdkp);
sp=sp+i*sp;
fm=zeros(nwin,K,nf,mdkp);
fm=fm+i*fm;
sv=zeros(nwin,nf,min([K,NCHAN]));
proj=zeros(N,K);
for j=1:nf 
    f0=f(j)/Fs;
    for k=1:K
      proj(:,k)=tapers(:,k).*exp(-f0*tvec');
    end
    tmp=data'*proj; % projected data
    [u,s,v]= svd(tmp,0); % svd 
    for mk=1:mdkp, 
      sp(n,:,j,mk)=u(:,mk)';
      fm(n,:,j,mk)=v(:,mk)';
    end  
    sv(n,j,:)=diag(s);
end;
