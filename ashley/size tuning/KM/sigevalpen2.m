function OF = sigevalpen2(c,x,y)

% double sigmoid fits Ae, ke=k1+k2, xe=x1, Ai, ki=k2, xi=x1+x2
% ke and xi defined so ex curve is steeper and in curve is centered higher
logfit2 = @(coefs,xdata) coefs(1)./(1+exp(-(coefs(2)+coefs(5))*(xdata-coefs(3)))) - coefs(4)./(1+exp(-coefs(5)*(xdata-(coefs(3)+coefs(6)))));

sse = norm(y - logfit2(c,x)).^2;

xrng = 0:0.1:max(x); % range from 0 to max size, step 0.1
f2 = diff(logfit2(c,xrng),2); % ~2nd derivative
pen = trapz(xrng(2:end-1),f2.^2);

lambda = 1e5;

OF = sse+lambda*pen;