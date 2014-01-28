function [best_dir, min_dir, R_best_dir, R_min_dir, sel_index] = find_best_stim(ydata)

% best_dir is the stim which elicited best response.
% min_dir is the stim which elicited least response.
% sel_index = 1 - R_min_dir / R_best_dir
%   Kenichi Ohki 04/01/08
%


nstim_per_run = length(ydata);
[R_best_dir, best_dir] = max(ydata);
[R_min_dir, min_dir] = min(ydata);
sel_index = 1-R_best_dir/R_min_dir;
