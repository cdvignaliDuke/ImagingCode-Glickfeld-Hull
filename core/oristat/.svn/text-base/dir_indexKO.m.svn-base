function [best_dir, null_dir, R_best_dir, R_null_dir, R_min_dir, DI] = dir_indexKO (ydata)

% classical direction index.
% DI = 1 - R(null_dir)/R(best_dir).
% null_dir is opposite direction to best_dir.
% best_dir is the direction which elicited best response.
%
%   Kenichi Ohki 09/20/04
%


nstim_per_run = length(ydata);
[R_best_dir, best_dir] = max(ydata);
R_min_dir = min(ydata);
null_dir = mod((best_dir + nstim_per_run/2 - 1), nstim_per_run)+1;
R_null_dir = ydata(null_dir);
DI = 1-R_null_dir/R_best_dir;
