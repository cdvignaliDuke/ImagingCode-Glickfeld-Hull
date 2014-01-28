function table = average_epoch (tc, epochs, nframes_per_run, run_ind)

% calculate average signal value during each epoch
%   (0deg grating, baseline, etc..) from single timecourse.
% 
% tc: timecourses of multiple ROIs. N x 1 or 1 x N.
% epochs: cell array (M x 1 or 1x M). Each cell represents one epoch 
%           (0 deg grating, 45 deg grating, baseline, etc...).
%           Each cell contains a series of numbers,
%           which indicate the frames for each epoch.
%           e.g. if epochs{1} represents 0 deg grating,
%           and 0 deg grating was displayed at 5-9 frames,
%           then epochs{1} = [5:9].
%           1-based (first frame is 1).
% nframes_per_run: number of frames per run
% run_ind: run indexes which you want to include in the table.
%           0-based (first run is 0).
% table: a table of (Number of runs) x(Number of epochs).
%           each cell contains an average value of timecourse
%           during the epoch.
%
%   Kenichi Ohki 09/16/04
%

Nframes = length(tc);

Nepochs = length(epochs);
Nruns = length(run_ind);

table = zeros(Nruns, Nepochs);

n=1;
for run=run_ind
    for i=1:Nepochs
        table(n, i)=sum(tc(run*nframes_per_run+epochs{i}))./length(epochs{i});
    end
    n=n+1;
end
