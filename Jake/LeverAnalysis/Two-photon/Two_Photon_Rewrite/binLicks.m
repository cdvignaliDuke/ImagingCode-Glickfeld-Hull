function [r_avg, sem_r] = binLicks(lick, ifi)
% Define number of columns to average
AVG_COLS = 3;

% Dimension over which to average
DIM = 2; % Columns

% Use filter to calculate the moving average across EVERY combination of columns
r_moving_avg = filter(ones(1,AVG_COLS)/AVG_COLS,1,lick,[],DIM);

% Grab only the column averages that were actually wanted
r_avg = r_moving_avg(:,AVG_COLS:AVG_COLS:end);
sem_r = std(r_avg,1)/sqrt(size(r_avg,1))/double(ifi)*1000;
r_avg = mean(r_avg,1)/double(ifi)*1000;
