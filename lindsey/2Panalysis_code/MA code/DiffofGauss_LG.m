%VB Gabor fitting, modified to fit Diff of Gauss
function h = DiffofGauss_LG(pars,Stims)
 a1 = pars(1);
 a2 = pars(2);
 b1 = pars(3);
 c1 = pars(4); 
 c2 = pars(5);
 R = pars(6);
 x = Stims;
 
 %h = a1*exp(-((x-b1)/c1).^2) - a2*exp(-((x-b1)/(c2)).^2);
 h = R + (a1.*c1.*erf((x-b1)./c1)) - (a2.*c2.*erf((x-b1)./c2));
 
%h = real(h);

