function [r,Hxy,Ht] = pinkrnd(pars,refreshRate)
% R = PINKRND(PARS,REFRESHRATE) Noise generator for PINKNOISE
% [R,H] = PINKRND(PARS,REFRESHRATE) also returns frequency spectrum
% PARS is a structure with the following fields
%  SEED random generator seed
%  NPIX number of pixels on a side
%  WIDTH stimulus size in degrees
%  LEN stimulus length in # refreshs
%  PAD # of blank frames before and after sequence
%  SF_C 
%  SF_0
%  TF_0
%  RANGE

% change log
%
% 09/07/23 vb fixed prctile bug
% 

if ~isfield(pars,'range')
    pars.range = 3;
end

rand('twister',pars.seed);
        
% spatial sampling frequency
Fxy = pars.npix/pars.width; % samples / deg

% spatial frequency grid
[ux,uy]=freqspace(pars.npix,'meshgrid'); % normalized frequency
vxy = Fxy/2*sqrt(ux.^2+uy.^2);
% temporal frequency vector
[ut,dummy]=freqspace(pars.len);       
vt = refreshRate/2*ut;

% amplitude spectrum, space
Hxy = (1./(vxy+pars.sf_c)).*abs(vxy<=pars.sf_0);
Hxy = repmat(Hxy,[1,1,pars.len]);

if pars.tf_0 > 0
    % amplitude spectrum, time
%     Ht = zeros(1,1,pars.len);
%     Ht(:) = abs(vt)<=pars.tf_0,[3 2 1];
%     Ht = repmat(Ht,[pars.npix,pars.npix,1]);

    Ht = permute(abs(vt)<=pars.tf_0,[3 1 2]);
    Ht = repmat(Ht,[pars.npix,pars.npix,1]);
else
    Ht = ones(1,1,pars.len);
    Ht = repmat(Ht,[pars.npix,pars.npix,1]);    
end

% space-time
H =Hxy.*Ht;

% Gaussian noise
x = randn(pars.npix,pars.npix,pars.len);

% pad with blank frames
x(:,:,1:pars.pad)=0;
x(:,:,end-pars.pad+1:end)=0;        

X = fftshift(fftn(x));       
Y = H.*X;
y = ifftn(fftshift(Y)); 

% VB fixed as of 09/07/21
% using PRCTILE is WRONG because it shifts DC around
% lims = prctile(y(:),[.1 99.9]);
% r = uint8((y*pars.c-lims(1))/(lims(2)-lims(1))*255);   
              
sigma = std(y(:));
lims = pars.range*[-sigma sigma];
r = uint8((y*pars.c-lims(1))/(lims(2)-lims(1))*255);

return;

%%
figure;
imagesq(uy)
imagesq(Hxy(:,:,1) .* ux>0 .* uy>0)

[th,rho]= cart2pol(ux,uy);

ang = pi/4;
wedge = exp(-((th+ang).^2)/2/.1^2) + exp(-((th+ang+pi).^2)/2/.1^2);
imagesq(wedge);

test = ifft2(fftshift(X(:,:,1).*Hxy(:,:,1) .* wedge));
imagesq(abs(test))

