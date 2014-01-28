function plotsig(C,sig,t,f)
% Function to plot C where it is higher than a threshold sig
% useful for plotting coherence
% Usage: plotsig(C,sig,t,f)
% Inputs:
% C: input array t x f. If vector, then as row vector
% sig: significance level
% t: t axis grid for plot
% f: f axis grid for plot.
if nargin < 4; error('Need all arguments'); end;
[T,F]=size(C);
if F==1; error('C needs to be a row vector'); end;
if T~=length(t) || F~=length(f);
    error('frequency and/or time axes are incompatible with data'); 
end;
if T==1;
    dim=max(T,F);
    C=C(:);
    mask=zeros(dim,1);
    indx=find(C>sig);
    mask(indx)=1;
    plot(f,mask.*C); 
    xlabel('f'); ylabel('|C|');
else
    mask=zeros(T,F);
    for n=1:length(t);
        for m=1:length(f);
           if C(n,m)>sig
              mask(n,m)=1;
           end;
        end;
    end;
    imagesc(t,f,(mask.*C)'); axis xy; colorbar
    xlabel('t'); ylabel('f');
end;