function rcimg = stackRevCor(stack,condlist,delays,Nframes_per_stim,run_ind)

% reverse correlate frames for each condition
% condlist: a list of condition numbers (1-based)
% delays: e.g. [-3:5], which means -3th frame to 5th frame from onset (=0). 
%   If you use Noff+Non stim structure, onset (frame 0) is the beginning of off-period (prestim off).  
% Nframes_per_stim: Noff + Non
% run_ind: 0-based, which run you want to include in the analysis
%
% rcimg: stack of ndelays x Nconds
%
%   June 24, 2009   Kenichi Ohki

if nargin < 5
    run_ind = 0;
end

Nconds=max(condlist); % number of different conditions
[xDim, yDim, tDim] = size(stack)
Nstim_per_run=length(condlist);
ndelays=length(delays);
Nframes_per_run = Nframes_per_stim * Nstim_per_run ;
Nrep=length(run_ind);
Nframes=Nframes_per_run * Nrep;



rcimg = zeros(xDim, yDim, ndelays * Nconds);

for cond = 1:Nconds
     cond
    condinds = find(condlist == cond);
    onset_frames = Nframes_per_stim * (condinds-1) + 1;
    for d = 1:ndelays
        fr = onset_frames + delays(d); 
        inds= []; 
        for i=1:Nrep
              inds = [inds;(run_ind(i)*Nframes_per_run + fr)];
        end
       inds
%        inds=(max(inds,1)); % do not sum over negative inds
        rcimgind = d + (cond-1)*ndelays;
        inds=mod(inds-1, Nframes)+1;
        rcimg(:,:, rcimgind) = sum(stack(:,:,inds),3)./length(inds);
 
    end
end

