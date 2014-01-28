function data=rmlinesc(data,params, p)
% removes significant sine waves from data (continuous data).
%
% Usage: data=rmlinesc(data,params, p)
%
%  Inputs:  
% Note that units of Fs, fpass have to be consistent.
%       data        (data in [N,C] i.e. time x channels/trials) - required.
%       params      structure containing parameters - params has the
%       following fields: tapers, Fs, fpass, pad
%	        tapers 	    (parameters for calculating tapers [NW,K]) - optional. Defaults to [3 5]
%	        Fs 	        (sampling frequency) -- optional. Defaults to 1.
%               fpass       (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%	        pad		    (padding factor for the FFT) - optional. Defaults to 0.  
%			      	    e.g. For N = 500, if PAD = 0, we pad the FFT 
%			      	    to 512 points; if PAD = 2, we pad the FFT
%			      	    to 2048 points, etc.
%	    p		    (P-value for F-test) - optional. Defaults to 0.05/N
%	    where N is data length. This corresponds to a false detect
%	    probability of approximately 0.05
%
%
%  Outputs: 
%       data        (data with significant lines removed)
%

[N,C]=size(data);
if nargin < 2 || isempty(params); params=[]; end;
[tapers,pad,Fs,fpass,err,trialave,params]=getparams(params);
if nargin<3 || isempty(p);p=0.05/N;end;
params.tapers=dpsschk(tapers,N,Fs); % calculate the tapers
[Fval,A,f,sig] = ftestc(data,params,p,'n');
fmax=findpeaks(Fval,sig);
for ch=1:C;
    fsig=f(fmax(ch).loc);
    Nf=length(fsig);
    fprintf('The significant lines for channel %d and the amplitudes are \n',ch);
    for nf=1:Nf;
        fprintf('%12.8f\n',fsig(nf));
        fprintf('%12.8f\n',real(A(fmax(ch).loc(nf),ch)));
        fprintf('%12.8f\n',imag(A(fmax(ch).loc(nf),ch))); 
        fprintf('\n');
    end;
    datasine(:,ch)=exp(i*2*pi*(0:N-1)'*fsig/Fs)*A(fmax(ch).loc,ch)+exp(-i*2*pi*(0:N-1)'*fsig/Fs)*conj(A(fmax(ch).loc,ch));
end;
% subplot(211); plot(data); hold on; plot(datasine,'r');
datan=data-datasine;
% subplot(212); plot(datan);
if nargout==0; 
   figure;subplot(211); plot(f,Fval); line(get(gca,'xlim'),[sig sig]);
   px1=pmtm(data(:,1));
   px2=pmtm(datan(:,1));
   subplot(212);plot(1:length(px1),10*log10(px1),1:length(px2),10*log10(px2));
end;
data=datan;   
