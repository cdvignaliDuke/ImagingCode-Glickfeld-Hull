function out = imups(im)
%IMUPS Upsample image by 2 using FFT.
% OUT = IMUPS(IM)

[m,n]=size(im);
mlarge=m*2;
nlarge=n*2;
out=zeros(mlarge,nlarge);
out(m+1-fix(m/2):m+1+fix((m-1)/2),n+1-fix(n/2):n+1+fix((n-1)/2)) = ...
    fftshift(fft2(im));

out = abs(ifft2(ifftshift(out)));

return;