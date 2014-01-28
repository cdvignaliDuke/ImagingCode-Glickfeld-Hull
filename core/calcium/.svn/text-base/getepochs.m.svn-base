function [stims,blanks,epochs] = getepochs (off, on, nstim, ntrials,pre,leica);
% [stims,blanks] = getepochs (off, on, nstim, ntrials);
%
% cell array {rep , stim}
%
% last stim is assigned to blank
%
% blank = entire blank period before each stim
%
% TODO: REWRITE
% make an epoch cell array from experiment parameters.
%   epochs for each stimulus (0 deg, 45 deg,...) and
%   an epoch for baseline are defined.
%   last cell contains frames for the baseline epoch.
% this is a recommendation. 
%   if you want to specify the epochs different from this,
%   please make a cell array for yourself.
%
% off: Number of off frames
% on: Number of on frames
% nstim: number of stimuli (usually gratings) per run
% pre: length (frames) of rest period before the beginnings of stimuli.
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

if nargin < 5
    pre = 1;
end
if nargin < 6
    leica = 0;
end

stims = cell(ntrials,nstim+1);
blanks = cell(ntrials,nstim);

trialdur = (off+on)*nstim;

nepochs = nstim + 1;

for itrial = 1:ntrials
    for istim=1:nstim
        % gattering off and on epoch for each stimulus
        blanks{itrial,istim} = (itrial-1)*trialdur+[(istim-1)*(off+on)+1 : istim*(off+on)-on];
        stims{itrial,istim} = (itrial-1)*trialdur+[(istim-1)*(off+on)+off+1 : istim*(off+on)];
        epochs{itrial,istim} = [blanks{itrial,istim} stims{itrial,istim}];
        
        if leica==1
            blanks{itrial,istim}=blanks{itrial,istim}-1;
            stims{itrial,istim}=stims{itrial,istim}-1;
            epochs{itrial,istim}=epochs{itrial,istim}-1;
        end
%        stims{itrial,istim} = stims{itrial,istim}(1:48);
        
        % gattering 
        if ~isempty(blanks{itrial,istim})
            stims{itrial,nepochs} = [stims{itrial,nepochs},blanks{itrial,istim}(end-pre+1:end)];
        end

    end
end

return;