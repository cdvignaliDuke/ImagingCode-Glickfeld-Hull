clear all

%% expt params
date = '110728';
mouse = 'AC41';
userun = [2];
nCond = 50;

nON = 19;
nOFF = 19;
nPlanes = 1;
pre_win = [11 19];
post_win = [20 38];

base = 'G:\users\lindsey\analysisLG\active mice';
outDir = fullfile(base, mouse,date);
%% resort expt
%load params and creat seq file
run(['PARAMS_' date '_' mouse]);
resort_seq_only

%load and concat stacks

stack = [];
for iRun = 1:length(userun);
    substack = readtiff(fullfile(base,mouse,date, [date '_' mouse '_run' num2str(userun(iRun)) '.tif']));
    stack = cat(3, stack, substack);
end
clear('substack');

%sort and average stacks
stack_sort
stack_sort_avg

%measure and sort by running
Running_LG