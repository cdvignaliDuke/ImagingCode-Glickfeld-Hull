function OF = sigevalpen1(c,x,y)

% single sigmoid fits Ae, ke=k1, xe=x1
logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))));

sse = norm(y - logfit1(c,x)).^2;

xrng = 0:0.1:max(x); % range from 0 to max size, step 0.1
f2 = diff(logfit1(c,xrng),2); % ~2nd derivative
pen = trapz(xrng(2:end-1),f2.^2);

lambda = 1e5;

OF = sse+lambda*pen;