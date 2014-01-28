function epochs = select_epochs (Noff, Non, nstim_per_run, rest_length);

% make an epoch cell array from experiment parameters.
%   epochs for each stimulus (0 deg, 45 deg,...) and
%   an epoch for baseline are defined.
%   last cell contains frames for the baseline epoch.
% this is a recommendation. 
%   if you want to specify the epochs different from this,
%   please make a cell array for yourself.
%
% Noff: Number of off frames
% Non: Number of on frames
% nstim_per_run: number of stimuli (usually gratings) per run
% rest_length: length (frames) of rest period before the beginnings of stimuli.
% epochs: cell array (N x 1). Each cell represents one epoch 
%           (0 deg grating, 45 deg grating, baseline, etc...).
%           Each cell contains a series of numbers,
%           which indicate the frames for each epoch.
%           e.g. if epochs{1} represents 0 deg grating,
%           and 0 deg grating was displayed at 5-9 frames,
%           then epochs{1} = [5:9].
%           1-based (first frame is 1).
%
%   Kenichi Ohki 09/16/04
%

Nepochs = nstim_per_run + 1;
epochs = cell(Nepochs,1);

for i=1:nstim_per_run
    epochs{i} = [(i-1)*(Noff+Non)+Noff : i*(Noff+Non)-1];
    epochs{Nepochs} = [epochs{Nepochs},...
            [(i-1)*(Noff+Non)+Noff-rest_length : (i-1)*(Noff+Non)+Noff-1]]; 
end



