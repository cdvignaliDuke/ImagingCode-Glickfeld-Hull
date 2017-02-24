%% Synthesize some data.
%
% play arround with noise_sigma to titrate the amount of "noise correlations"
%

clear all; close all; clc

N_time = 2; % seconds
N_trials = 1e3;
dt = 0.001;
tt = 0:dt:N_time;

noise_sigma = 1;
noise_mean = 0;


sin_freq = 5; % hz


% make a signal template to add to each population. each row is a trial.
% This will be the thing that creates "noise" correlations between cells.
% The correlated noise will be a sin wave. The phase will be randomized
% from trial to trial so that we can distinguish true trial-by-trial
% correlations from other types of correlations
phase = unifrnd(0, 2*pi, N_trials, 1);
tt_mtx = repmat(0:dt:N_time, N_trials, 1);
innerpart = 2.*pi.*tt_mtx.*sin_freq;
innerpart = bsxfun(@plus, innerpart, phase);
sin_waves = sin(innerpart); % plot this to verify that each row is a sin wave


% make matrix for "cell one". each row is a trial. cell_one responses go
% up, but have some correlated noise with cell_two
m = 2;
b = 4;
cell_one = (m.*tt_mtx + b) + normrnd(noise_mean, noise_sigma, size(tt_mtx)) + sin_waves;


% make matrix for "cell two". cell_two responses go down, but have some
% correlated noise with cell_one
m = -2;
b = 0;
cell_two = (m.*tt_mtx + b) + normrnd(noise_mean, noise_sigma, size(tt_mtx)) + sin_waves;



%% calculate correlations via CORRCOEF

% if you give 'corrcoeff' matricies as inputs, the funtion un-wraps the
% matrix to form column vectors. Pre-process the original data matricies
% [N_trials x N_time]. Subtract off the average time vector, then
% transpose so that un-wrapping works correctly.
%
% it's important to subtract off the average "time-series" and not the
% "averages across time" for each trial
%
cell_1_preprocess = transpose(bsxfun(@minus, cell_one, mean(cell_one, 1)));
cell_2_preprocess = transpose(bsxfun(@minus, cell_two, mean(cell_two, 1)));
corr_mtx =  corrcoef(cell_1_preprocess, cell_2_preprocess);
corr_from_corrcoef = corr_mtx(1,2)

%% calculate correlation via dot products

% always important to subtract off the average time-series, and transpose
% so time goes down columns.
cell_1_preprocess = transpose(bsxfun(@minus, cell_one, mean(cell_one, 1)));
cell_2_preprocess = transpose(bsxfun(@minus, cell_two, mean(cell_two, 1)));

cell_1_preprocess = cell_1_preprocess(:);
cell_2_preprocess = cell_2_preprocess(:);

% notice the transpose, and the * (not .*), and the ./ (not /)
corr_from_dotprod = (cell_1_preprocess' * cell_2_preprocess) ./ (norm(cell_1_preprocess) .* norm(cell_2_preprocess))


%% comapare to CORR
% use CORR to calculate all pair-wise correlations among trials. This
% process is expensive because the off-diagonal elements of the correlation
% matrix should not covary on a trial-by-trial basis (I randomized the
% phase). So calling CORR calculates gazillions more things than we need
% to, and takes lots longer.
%
% But we can calculate the distribution of correlations on a trial-by-trial
% basis, and verify that when we break the correspondance between trials
% the correlations go away (phases were randomized). The average
% correlation for shuffled trials should be zero

% always important to subtract off the average time-series, and transpose
% so time goes down columns.
cell_1_preprocess = transpose(bsxfun(@minus, cell_one, mean(cell_one, 1)));
cell_2_preprocess = transpose(bsxfun(@minus, cell_two, mean(cell_two, 1)));

all_corrs = corr(cell_1_preprocess, cell_2_preprocess);
avg_corr_from_full_corr_mtx = mean(diag(all_corrs))

figure
subplot(1,2,1)
histogram(diag(all_corrs))
xlabel('correlation')
title('trial-by-trial correlations')
subplot(1,2,2)
off_diag = triu(all_corrs, 1);
off_diag = off_diag(:);
off_diag(off_diag==0) = [];
histogram(off_diag)
xlabel('correlation')
title('correlations when trials shuffled')


