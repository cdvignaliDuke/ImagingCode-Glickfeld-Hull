function [xycov,lags] = xcovu(x,y,varargin)
% [XYCOV,LAGS] = XCOVU(X,Y,VARARGIN)

[xycov,lags] = xcov(x,y,varargin{:},'unbiased');

xycov = xycov /std(x)/std(y);

return

%% note that xcov normalizes for autocorrelation

x1 = randn(1,100000);
x2 = randn(1,100000);

y = filter(gausswin(8),1,x1)+.1*x2;

clf
subplot(2,1,1);
plot(xcovu(y,x2,64),'-k.'); 
subplot(2,1,2);
plot(xcovu(y,x1,64),'-k.'); 