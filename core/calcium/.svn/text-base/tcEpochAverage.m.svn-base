function [av,sd] = tcEpochAverage (timeCourses, epochs, ind)
%TCEPOCHAVERAGE Calculate average signal over epochs.
% av = TCEPOCHAVERAGE (TIMECOURSES, EPOCHS) where
%  TIMECOURSES is time courses of multiple ROIs. N x 1 or 1 x N.
%  EPOCHS is cell array (M x 1 or 1x M). Each cell represents one epoch 
%           (0 deg grating, 45 deg grating, baseline, etc...).
%           Each cell contains a series of numbers,
%           which indicate the frames for each epoch.
%           e.g. if epochs{1} represents 0 deg grating,
%           and 0 deg grating was displayed at 5-9 frames,
%           then epochs{1} = [5:9].
%           1-based (first frame is 1).
%
% av is av of (Number of runs) x(Number of epochs).
%           each cell contains an average value of timecourse
%           during the epoch.
%
%   Kenichi Ohki 09/16/04
%

if nargin < 3
    ind = [];
end

[nsamples,ncells] = size(timeCourses);

[ntrials,nepochs]=size(epochs);

av = cell(ncells,1);
sd = cell(ncells,1);

for itrial = 1:ntrials
    for istim = 1:nepochs
        for icell = 1:ncells
            sel = epochs{itrial,istim};
            if ~isempty(ind)
                sel = sel(ind);
            end
            av{icell}(itrial,istim) = double(mean(timeCourses(sel,icell)));
            sd{icell}(itrial,istim) = double(std(timeCourses(sel,icell)));
        end
    end

end

return;
