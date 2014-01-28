function [r,coefs] = waveletrnd3(pars)
%WAVELETRND3 Wavelet noise generator for WAVELETNOISE3D.
% R = WAVELETRND3(PARS)
%
% This is the stimulus used for the paper entitled:
% Robust characterization of receptive fields in mouse visual cortex

if ~isfield(pars,'range')
    pars.range = 3;
end

if ~isfield(pars,'Js')
    pars.Js = [];
end

total = 0;

[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;

J = log2(pars.npix)-1; % number of scales (excl lowest freq)
rand('twister',pars.seed);

% generate 4x the area and crop to avoid wrapping
coefs = dualtree3dcoefs([pars.npix*2,pars.npix*2,pars.len],J,pars.p0,pars.k,[pars.Jmin pars.Jmax],pars.Js);

inverse  = idualtree3D(coefs, J, Fsf, sf);

y = inverse(1:pars.npix,1:pars.npix,:); % crop

% VB fixed as of 09/07/21
% using PRCTILE is WRONG because shifts DC around
% lims = prctile(y(:),[.5 99.5]);
% r = uint8((y*pars.c-lims(1))/(lims(2)-lims(1))*255);   
 
sigma = std(y(:));
lims = pars.range*[-sigma sigma];
r = uint8((y*pars.c-lims(1))/(lims(2)-lims(1))*255);

return;
