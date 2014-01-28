function [w,J] = dualtree2dcoefs(dims,J,p0,k,Jrange)
%DUALTREE3DCOEFS
% [W,J] = DUALTREE2DCOEFS(DIMS,J,P0,K,JMIN,SEED)
% where DIMS = [NY,NX]
%       J is number of spatial scales
%       P0 is probability of nonzero coefficients
%       K is exponent of power-law decay of amplitude of coefficients
%       JMIN is the fines spatial scale
%       SEED is seed of random generator, rand('twister',seed)

if nargin < 1
    dims = [32,32];
end

% set seed
ny = dims(1);
nx = dims(2);

if nargin < 2
    J  = log2(ny)-1; % number of scales;
end

if nargin < 3
    p0 = .99;
end

if nargin < 4
    k = 2;
end

if nargin < 5
    Jrange = [2 J+1];
end

% probabilities 
probs = 1-(1-p0)*2.^[0:J-1];

% (1-p0) % probability of coefficients being non zero
% sum(ncoefs)*(1-p0)./ncoefs
% amplitude proportional to spatial period to the power k
amplitudes = ([2.^[0:J]].^k); 

w={};

% For j = 1..J, k = 1..2, d = 1..3, w{j}{k}{d} 
% are the wavelet coefficients produced at scale j 
% and orientation (k,d).


% high frequencies
for indJ = 1:J % from fine to coarse
    for indI = 1:2
        for indD = 1:3
            siz = [ny,nx]/2^(indJ);
            w{indJ}{indI}{indD} = zeros(siz);
            
            if indJ >= Jrange(1) && indJ <= Jrange(2)
                
                % decide which coefficient are nonzero
                rolls = find(rand(siz)>probs(indJ));
                
                % decide which coeffcient are negative
                heads = rand(size(rolls))<.5;
                
                pos = find(heads);
                neg = find(~heads);
                
                w{indJ}{indI}{indD}(rolls(pos)) = amplitudes(indJ);
                w{indJ}{indI}{indD}(rolls(neg)) = -amplitudes(indJ);

            end
        end
    end
end
    
siz = [ny,nx]/2^(J);

% low frequencies
for indI = 1:2
    w{J+1}{indI} = zeros(siz);

    if J+1<= Jrange(2)
        rolls = find(rand(siz)>probs(J));
        heads = rand(size(rolls))<.5;
        pos = find(heads);
        neg = find(~heads);
        w{J+1}{indI}(rolls(pos)) = amplitudes(indJ);
        w{J+1}{indI}(rolls(neg))= -amplitudes(indJ);
    end
end

return;

%%
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
[w,J] = dualtree2dcoefs;
y = idualtree2D(w, J, Fsf, sf);
imagesc(y);axis square; colormap gray;
