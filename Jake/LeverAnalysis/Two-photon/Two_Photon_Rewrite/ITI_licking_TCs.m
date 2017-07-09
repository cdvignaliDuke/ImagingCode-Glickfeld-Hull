%script for using ITI lick bouts to ID licking related neurons and ID them
%will us dfoverf_TC

%THINGS TO DO BEFORE THIS CODE IS USEABLE
%1) figure out counter alignment issue
%2) apply mask to all useable frames for animals limiited to the first 70k
%3) finish the actual code


% Window shoulld start 5s after reward delivery and end 4s before cue onset
%timepoints needed: reward delivery and cue onset
%want to plot 1s before lick bout onset and 4s after onset. 

% need to make sure OR and UR trials do not fuck this up
%trial_outcome variables zeroed to .....
%for CRP trials lever.press is the time of cue onset
%unequal # of frames in frame_info, lick_data and df/f TC

clear;
file_info_CRP;
dir_base = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\';

for ii = [2, 5, 9, 22, 23, 31]; %1:length(dates)
    %assign an output directory for the results
    out_dir = [dir_base, 'ITI licking\', dates{ii}, '_', mouseID{ii}];
    
    %check to see if the file exists
    for rID = 1:3;
        data_dir = [dir_base, dates{ii}, '_', runID{rID}, '_', mouseID{ii}, '\'];
        TC_fn = dir([data_dir, '_dFOverF_TC.mat']);
        if ~isempty(TC_fn);
            break
        end
    end
    
    %load TCs  bx data and licking data
    load([data_dir, TC_fn.name]);
    load([data_dir, 'parse_behavior'], 'lick_data', 'lever');
    
end