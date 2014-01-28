function [w,J] = waveletcoefs4(dims,J,p0,Jc)
%WAVELETCOEFS4
% [W,J] = WAVELETCOEFS4(DIMS,J,P0,K,JRANGE,JS)
% where DIMS = [NY,NX,NT]
%       J is number of spatial scales (excluding lowest freq)
%       SEED is seed of random generator, rand('twister',seed)
%       P0 is probability of nonzero coefficients
%       K is exponent of power-law decay of amplitude of coefficients
%       JRANGE = [JMIN,JMAX] are the finest and lowest spatial scales
%       JS is scale to output, default = [] (all scales). this is a hack to
%       look at one scale at time.

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

% probability as function of scale
probs = 1-(1-p0)*(2.^[0:J-1]);

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
            
            if indJ == Jc
                  
                % decide which coefficient are nonzero
                rolls = find(rand(siz)>probs(indJ));
                
                % decide which coeffcient are negative
                heads = rand(size(rolls))<.5;
                
                pos = find(heads);
                neg = find(~heads);
                                
                w{indJ}{indI}{indD}(rolls(pos)) = 1;
                w{indJ}{indI}{indD}(rolls(neg)) = -1;

            end;
        end;
    end;
end;
    
siz = [ny,nx,nt]/2^(J);

% low frequencies
for indI = 1:4
    w{J+1}{indI} = zeros(siz);

    if J+1 == Jc 
        rolls = find(rand(siz)>probs(J));
        heads = rand(size(rolls))<.5;
        pos = find(heads);
        neg = find(~heads);
        
        w{J+1}{indI}(rolls(pos)) = 1;
        w{J+1}{indI}(rolls(neg))= -1;
    end
end

return;

%%
[Faf, Fsf] = FSfarras;
[af, sf] = dualfilt1;
y = idualtree3D(w, J, Fsf, sf);

%% 

