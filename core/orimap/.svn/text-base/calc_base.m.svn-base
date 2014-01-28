function [base, base_sm] = calc_base (ave, base_inds, nframes_per_stim, nstim_per_run, sp_filter)

% calculate a baseline image
% base_inds specify the baseline frames in one trial. e.g. [1:4].
%   one trial starts from frame 1 and end at frame 'nframes_per_stim'.
% nframes_per_stim: number of frames in one trial (=stimulation).
% nstim_per_run: number of stims in one run.
% sp_filter: spatial low-pass filter
%   it this parameter is omitted, base_sm=0.
% base: baseline image without smoothing
% base_sm: baseline image with smoothing

dim=size(ave);
base=zeros(dim(1),dim(2));

for n=1:nstim_per_run   
    for i=base_inds
        base=base+ave(:,:,(n-1)*nframes_per_stim+i);
    end
end
base=base./(length(base_inds)*nstim_per_run);

if (nargin >=5)
    base_sm=filter2(sp_filter, base);
else
    base_sm=0;
end
