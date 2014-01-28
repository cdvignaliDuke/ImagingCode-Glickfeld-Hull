function [params]  = imFindCellsParamsSet(varargin)
%IMFINDCELLSPARAMSSET Creates / alters parameters structure for IMFINDCELLS
%   [params]  = imFindCellsParamsSet(varargin)
% s  structure (or array of structures) with filenames etc.
% nf:  index of array of structures.  Optional.  If s is a single structure than function is called with single argument.

% default parameters
params.min_area = 30;
params.win_size = 45; % not in use
params.noise_radius = 1; 
params.junk_radius = 2;
params.ratio_th = 0.23;
params.min_center = 25;
params.con_ratio = 0.65;
params.breakup_factor = 3;
params.clear_border = 0;
params.fast_breakup = 0;
params.do_manual = 0;
params.show_flag = 1;

index = 1;
while index < length(varargin)-1
    params.(varargin{index})=varargin{index+1};
    index = index + 2;
end

return;
