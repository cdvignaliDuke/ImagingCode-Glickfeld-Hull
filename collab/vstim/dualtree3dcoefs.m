function [w,J] = dualtree3dcoefs(dims,J,p0,k,Jrange,Js)
%DUALTREE3DCOEFS
% [W,J] = DUALTREE3DCOEFS(DIMS,J,P0,K,JRANGE,JS)
% where DIMS = [NY,NX,NT]
%       J is number of spatial scales (excluding lowest freq)
%       SEED is seed of random generator, rand('twister',seed)
%       P0 is probability of nonzero coefficients
%       K is exponent of power-law decay of amplitude of coefficients
%       JRANGE = [JMIN,JMAX] are the finest and lowest spatial scales
%       JS is scale to output, default = [] (all scales). this is a hack to
%       look at one scale at time.
%
% This is the stimulus used for the paper entitled:
% Robust characterization of receptive fields in mouse visual cortex


if nargin < 1
    dims = [32,32,1024];
end

% set seed
ny = dims(1);
nx = dims(2);
nt = dims(3);

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

%rand('twister',seed);

% high frequencies
for indJ = 1:J % from fine to coarse
    for indI = 1:4
        for indD = 1:7
            siz = [ny,nx,nt]/2^(indJ);
            w{indJ}{indI}{indD} = zeros(siz);
            
            if indJ >= Jrange(1) && indJ <= Jrange(2)
                
                % decide which coefficient are nonzero
                rolls = find(rand(siz)>probs(indJ));
                
                % decide which coeffcient are negative
                heads = rand(size(rolls))<.5;
                
                pos = find(heads);
                neg = find(~heads);
                
                % to output one spatial scale at a time
                if ~isempty(Js)
                    select = double(indJ == Js);
                else
                    select = 1;
                end
                
                w{indJ}{indI}{indD}(rolls(pos)) = amplitudes(indJ) * select;
                w{indJ}{indI}{indD}(rolls(neg)) = -amplitudes(indJ) * select;

            end
        end
    end
end
    
siz = [ny,nx,nt]/2^(J);

% low frequencies
for indI = 1:4
    w{J+1}{indI} = zeros(siz);

    if J+1<= Jrange(2)
        rolls = find(rand(siz)>probs(J));
        heads = rand(size(rolls))<.5;
        pos = find(heads);
        neg = find(~heads);
        
        % to output one spatial scale at a time
        if ~isempty(Js)
            select = double(J+1 == Js);
        else
            select = 1;
        end
                
        w{J+1}{indI}(rolls(pos)) = amplitudes(indJ) * select;
        w{J+1}{indI}(rolls(neg))= -amplitudes(indJ) * select;
    end
end

return;

%%
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
y = idualtree3D(w, J, Fsf, sf);

%% 

