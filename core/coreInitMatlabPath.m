function coreInitMatlabPath(svnroot,ijroot)
%COREINITMATLABPATH  Initializes MATLAB for Reid Lab Core
% COREINITMATLABPATH(SVNROOT,IJROOT)

svncore = fullfile(svnroot,'core');

addpath(svncore);
addpath(fullfile(svncore, 'tools'));
addpath(fullfile(svncore, 'calcium'));
addpath(fullfile(svncore, 'fastrig'));
addpath(fullfile(svncore, 'multicore'));
addpath(fullfile(svncore, 'reg'));
addpath(fullfile(svncore, 'calciumsandbox'));
addpath(fullfile(svncore, 'oristat'));
addpath(fullfile(svncore, 'orimap'));


% allow assertions for old versions of matlab
if ~(exist('assert') == 5)  % builtin function
    addpath(fullfile(svncore, 'tools', 'assert-old'));
end

%disp('COREINITMATLABPATH: Added a bunch of dirs to the matlab path.');

return;
