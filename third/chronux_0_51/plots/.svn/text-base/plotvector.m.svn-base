function plotvector(X,f,plt,Xerr)
% Function to plot a frequency dependent vector X. If error bars are specified in Xerr,
% it also plots them. Xerr can either contain upper and lower confidence intervals 
% on X, or simply a theoretical confidence level (for the coherence). Used
% to plot the spectrum and coherency.
% Usage: plotvector(X,f,plt,Xerr)
% Inputs:
% X: input vector as a function of frequency (f), see third argument
% f: f axis grid for plot. Default. [1:length(X)]
% plt: 'l' for log, 'n' for no log.
% Xerr: lower and upper confidence intervals for X1: lower/upper x f. Or
%       simply a single number specifying an f-independent confidence
%       level.
if nargin < 1; error('Need data'); end;
N=length(X); 
if nargin < 2;
    f=1:N;
end;
if length(f)~=N; error('frequencies and data have incompatible lengths'); end;
if nargin < 3;
    plt='l';
end;
if strcmp(plt,'l');
    X=10*log10(X);
    if nargin ==4; Xerr=10*log10(Xerr); end;
end;

if nargin < 4;
    plot(f,X);
else
    if length(Xerr)==1;
       plot(f,X); 
       line(get(gca,'xlim'),[Xerr,Xerr],'Color','r');
    else
       plot(f,X); 
       hold on; plot(f,Xerr(1,:),'g'); plot(f,Xerr(2,:),'g'); 
    end
end
xlabel('f');
if strcmp(plt,'l'); ylabel('10*log10(X)'); else ylabel('X'); end;


