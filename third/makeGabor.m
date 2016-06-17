function gabor = makeGabor(imSize, lambda, theta, phase, contrast);
% make sinusoidal gabor
% imSize: dimensions in pixels nxn
% lambda: pixels/cycles
% theta: orientation
% phase: 0 to 1;
% contrast: 0 to 1

sigma = 10;                             % gaussian standard deviation in pixels
trim = .01 ./ contrast+eps;                 % trim off gaussian values smaller than this

X = 1:imSize;                           % X is a vector from 1 to imageSize
X0 = (X / imSize) - .5;                 % rescale X -> -.5 to .5
%figure;                                 % make new figure window
%plot(X0);                               % plot ramp

sinX = sin(X0 * 2*pi);                  % convert to radians and do sine
%plot(sinX);                             % plot a 1D sinewave!

freq = imSize/lambda;                    % compute frequency from wavelength
Xf = X0 * freq * 2*pi;                  % convert X to radians: 0 -> ( 2*pi * frequency)
sinX = sin(Xf) ;                        % make new sinewave
% plot(sinX, 'r-');                       % plot in red
phaseRad = (phase * 2* pi);             % convert to radians: 0 -> 2*pi
sinX = sin( Xf + phaseRad) ;            % make phase-shifted sinewave
% hold on;                                % superimpose next plot on last
% plot(sinX, 'g-');                       % plot in green
% hold off;                               % next plot overwrites this one

[Xm Ym] = meshgrid(X0, X0);             % 2D matrices
%imagesc( [ Xm Ym ] );                   % display Xm and Ym
%colorbar; axis off                      % add colour bar to see values
Xf = Xm * freq * 2*pi;
grating = sin( Xf + phaseRad);          % make 2D sinewave
% imagesc( grating, [-1 1] );             % display
% colormap gray(256);                     % use gray colormap (0: black, 1: white)
% axis off; axis image;                   % use gray colormap

thetaRad = (theta / 360) * 2*pi;        % convert theta (orientation) to radians
Xt = Xm * cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym * sin(thetaRad);                % compute proportion of Ym for given orientation
XYt = [ Xt + Yt ];                      % sum X and Y components
XYf = XYt * freq * 2*pi;                % convert to radians and scale by frequency
grating = sin( XYf + phaseRad);                   % make 2D sinewave
% imagesc( grating, [-1 1] );                     % display
% axis off; axis image;    % use gray colormap

s = sigma / imSize;                     % gaussian width as fraction of imageSize
Xg = exp( -( ( (X0.^2) ) ./ (2* s^2) ));% formula for 1D gaussian
Xg = normpdf(X0, 0, (20/imSize)); Xg = Xg/max(Xg);  % alternative using normalized probability function (stats toolbox)
% plot(Xg)

gauss = exp( -(((Xm.^2)+(Ym.^2)) ./ (2* s^2)) ); % formula for 2D gaussian
% imagesc( gauss, [-1 1] );                        % display
% axis off; axis image;     % use gray colormap
gauss(gauss < trim) = 0;                 % trim around edges (for 8-bit colour displays)
gabor = grating .* gauss .* contrast;                % use .* dot-product
% imagesc( gabor, [-1 1] );                        % display
% axis off; axis image;                    % use gray colormap
% axis image; axis off; colormap gray(256);
% set(gca,'pos', [0 0 1 1]);               % display nicely without borders
% set(gcf, 'menu', 'none', 'Color',[.5 .5 .5]); % without background