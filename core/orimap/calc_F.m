function dir_F = calc_F (ave, stim_inds, nframes_per_stim, nstim_per_run)

% calculate images to each stimulation
% stim_inds specify the activation frames in one trial. e.g. [5:9].
%   one trial starts from frame 1 and end at frame 'nframes_per_stim'.
% nframes_per_stim: number of frames in one trial (=stimulation).
% nstim_per_run: number of stims in one run.

dim=size(ave);
dir_F=zeros(dim(1),dim(2),nstim_per_run);

for n=1:nstim_per_run
    for i=stim_inds
        dir_F(:,:,n)=dir_F(:,:,n)+ave(:,:,(n-1)*nframes_per_stim+i);
    end
end
dir_F=dir_F./length(stim_inds);