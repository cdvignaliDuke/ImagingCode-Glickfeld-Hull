function out=imHueRotateKO (in, dhue)

% rotate hue of an image in KO space
% dhue: 0-1
%
% Kenichi Ohki 2008.11.4

Nx=size(in,1);
Ny=size(in,2);

temp=reshape(in,Nx*Ny,3);
hsv=rgb2hsv(temp);
h=hsv(:,1);
hKO=hueHSV2KO(h);
hKOrotated=mod(hKO+dhue,1);
hrotated=hueKO2HSV(hKOrotated);
hsv(:,1)=hrotated;
rgb=hsv2rgb(hsv);
out=reshape(rgb,Nx,Ny,3);







