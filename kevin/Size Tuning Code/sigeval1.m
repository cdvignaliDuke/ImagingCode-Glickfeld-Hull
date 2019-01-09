function OF = sigeval1(c,x,y)

% single sigmoid fits Ae, ke=k1, xe=x1
logfit1 = @(coefs,xdata) coefs(1)./(1+exp(-coefs(2)*(xdata-coefs(3))));

sse = norm(y - logfit1(c,x)).^2;

OF = sse;