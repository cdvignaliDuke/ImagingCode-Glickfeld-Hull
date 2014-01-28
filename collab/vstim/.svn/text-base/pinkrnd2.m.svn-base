function [r,Hxy,Ht] = pinkrnd2(pars,refreshRate)
% R = PINKRND2(PARS,REFRESHRATE) Noise with 
% spatial spectrum ~1/(f + fc)*( f < f0)
% and temporal spectrum 
% [R,H] = PINKRND2(PARS,REFRESHRATE)
% PARS is a structure with the following fields
%  LEN stimulus length in # refreshs
%  PAD # of blank frames before and after sequence
%  NPIX number of pixels on a side
%  C contrast
%  SF_C 
%  SF_0
%  TF_0
%  SEED random generator seed
%  WIDTH stimulus size in degrees
%  RANGE
%  EQ if 1 performs histeq on stimulus stack

if ~isfield(pars,'thresh');
    pars.thresh = 0;
end

if ~isfield(pars,'range')
    pars.range = 3;
end

if ~isfield(pars,'eq')
    pars.eq = 0;
end

rand('twister',pars.seed);
        
% spatial sampling frequency
Fxy = pars.npix/pars.width; % samples / deg

% spatial frequency grid
[ux,uy]=freqspace(pars.npix,'meshgrid'); % normalized frequency

vrho = Fxy/2*sqrt(ux.^2+uy.^2); % cycles / deg
% temporal frequency vector
[ut,dummy]=freqspace(pars.len);       
vt = refreshRate/2*ut; % Hz

%% design filter

% sigma chosen so that window does not blur
% desired filter frequency characteristics while
% reducing ripples in space
sigma_pix = pars.npix/ 12;

w = fspecial('gaussian',pars.npix,sigma_pix);
w = w ./ max(w(:));  % Make the maximum window value be 1.

% ideal frequency response
Ixy = (1./(vrho+pars.sf_c)).*abs(vrho<=pars.sf_0);
%Ixy = abs(vrho<=pars.sf_0);

hxy = fwind2(Ixy,w);
Hxy = fft2(hxy);

%% illustrates how ripples eliminated
% 
% figure;
% subplot(321);
% imagesq(fftshift(abs(fft2(Ixy))));
% title('Ideal filter (space)');
% subplot(322);
% imagesq(Ixy)
% title('Ideal filter (frequency)');
% subplot(324);
% imagesq(fftshift(abs(Hxy)));
% title('Windowed filter (frequency)')
% subplot(323);
% imagesq(hxy);
% title('Windowed filter (space)');
% subplot(325);
% disc = abs(vrho<=.25);
% imagesq(disc,[-1 1]);
% title('Example image');
% subplot(326);
% imagesq(fftshift(ifft2(fft2(disc).*Hxy)),[-1 1]);
% axis tight;
% title('Filtered image');
% supertitle(sprintf('Spatial filter design \\sigma = %2.1f ',sigma_pix),.98);

%%
%x = randn(pars.npix,pars.npix,pars.len);
x = rand(pars.npix,pars.npix,pars.len) - .5;

% pad with blank frames
x(:,:,1:pars.pad)=0;
x(:,:,end-pars.pad+1:end)=0;        

X = fftn(x);
Y = bsxfun(@times,X,Hxy);
%y = fftshift(ifft2(Y));

% cutoff = refreshRate/sigma;
% alpha = 1/sigma;
% alpha = cutoff/refreshRate; % cut-off normalized frequency

% how to design gaussian filter
% alpha = pars.tf_0 / (refreshRate/2);
% n = 64
% stem(freqspace(n,'whole')*refreshRate/2,(abs(fft(gausswin(n,n*alpha)))))

% gaussian filter with cut-off tf_0
% alpha = pars.tf_0 / (refreshRate/2);
% Ht = permute(fft(gausswin(pars.len,pars.len*alpha)),[3 2 1]);

% pars.tf_0 = 4
b = fir1(pars.len-1,pars.tf_0/(refreshRate/2));

Ht = permute(fft(b),[3 1 2]);

% plot(vt,fftshift(abs(squeeze(Ht)))/max(fftshift(abs(squeeze(Ht)))))
% hold on;
% grid

Y2 = bsxfun(@times,Y,Ht);

% fftshift needed here because of some quirks I don't
% understand about the way spatial filter is designed
y = fftshift(ifftn(Y2));

sigma = std(y(:));

lims = pars.range*[-sigma sigma];

% interval [0,1]
s = (y*pars.c-lims(1))/(lims(2)-lims(1));

if pars.thresh
    r = uint8((s>=.5)*255);
else
    r = uint8(s*255);
end

if pars.eq
    % make one 2d image
    r = reshape(r,[pars.npix*pars.npix,pars.len]);
    r = histeq(r);
    r = reshape(r,[pars.npix,pars.npix,pars.len]);
end

%%
return;
