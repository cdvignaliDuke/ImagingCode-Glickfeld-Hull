%VB Gabor fitting, modified to fit 2D elliptical
function h = Gauss2D_ellipseMA(pars,sftf)
%H = 2DGauss_ellipseMA(PARS,XX,YY)
%
% %assume no skew.. 
% A = pars(1);
% sigma_SF = pars(2); % in cycles/deg
% sigma_TF = pars(3); % in Hz
% sf0 = pars(4); %center sf
% tf0 = pars(5); %center tf
%%% sf =SF; % cycles / deg
%%% tf = TF; % Hz

 A = pars(1);
 sigma_SF = pars(2); % in log2(cycles/deg)
 sigma_TF = pars(3); % in log2(Hz)
 log2sf = sftf(:,1); % log2(cycles / deg)
 log2tf = sftf(:,2); % log2(Hz)
 log2sf0 = pars(4); %center log2(sf)
 log2tf0 = pars(5); %center log2(tf)
 xi = pars(6);
 log2tfpsf = xi*(log2sf - log2sf0) + log2tf0; 
 
 h = A*exp(-1*((log2sf - log2sf0).^2)./(2*sigma_SF^2)) .* exp(-1*((log2tf -  log2tfpsf).^2)./(2*sigma_TF^2));
 

%h = real(h);

