function [r,coefs] = waveletrnd4(pars)
%WAVELETRND4 Wavelet noise generator for WAVELETNOISE3D.
% R = WAVELETRND4(PARS)
%
% This is the stimulus used for the paper entitled:
% Robust characterization of receptive fields in mouse visual cortex


if ~isfield(pars,'range')
    pars.range = 3;
end;

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;

% number of spatial scales 
J = log2(pars.npix)-1; 
rand('twister',pars.seed);

%% 
y = zeros(pars.npix,pars.npix,pars.len);
nIter = pars.Jmax - pars.Jmin + 1;

% generate one scale at a time, downsample in space to keep temporal
% spectrum constant,iterate from large to small scales

for iter = 1:nIter
    fprintf('Scale %i/%i\n',iter,nIter);
    os = 2^(iter-1);
    npix = pars.npix*os*2;
    J =  log2(npix)-1;
    w = waveletcoefs4([npix,npix,pars.len], J, pars.p0, pars.Jmax);
    W  = idualtree3D(w, J, Fsf, sf);
    W2 = imresize(W,1/os); % downsample in space

    y = y + double(W2(1:pars.npix,1:pars.npix,:)); % crop
end;
 
sigma = std(y(:));
lims = pars.range*[-sigma sigma];
r = uint8((y*pars.c-lims(1))/(lims(2)-lims(1))*255);

%%
return;

