function reconstructTestAll
%
% Call reconstructTestEach on several different series.
% This is meant to be a unit test for ijdad2tiffseq: a quick way to check
% for successful reconstruction
%
% 080619 histed: v1
%
% $Id$

dbstop if error;

%% version
% $$$ addpath('i:/users/histed/matlab/reidlab-old-versions/fastrig-2008apr10', ...
% $$$         '-begin');

%% consts
dadDataRoot = '//zmey/storlab/data/vincent';
outRoot = 'i:/users/histed/test-output';
exptList = { 'mouse080327', 'estimB'; ...
             'mouse080327', 'estim2B'; ...             
             'mouse080327', 'spont60mWA'; ...
             'mouse080612', 'straight150mW';
             };


%% iterate over the list
nExpts = size(exptList,1);
for iE=1:nExpts
    tExptName = exptList{iE,1};
    tSeriesName = exptList{iE,2};    
    reconstructTestEach(dadDataRoot, outRoot, ...
                        tExptName, tSeriesName);
end
