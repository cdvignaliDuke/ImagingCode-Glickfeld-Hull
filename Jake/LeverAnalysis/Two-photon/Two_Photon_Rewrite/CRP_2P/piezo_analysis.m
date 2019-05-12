%Function for loading and analyzing piezo data
%will need to integrate into the main analysis pipeline
%will need to insert into getTC_events_CRP somewhere. 
%should make another function which is a new version of
%trigger_movie_by_event_2P which extracts piezo traces for each trial

%this current one should merely retreive the full piezo trace and align it
%to the imaging frames

function [mov_data, mov_frame_ind] = get_piezo_data(move_data_dir, mov_file);

mov_data_temp = fopen([move_data_dir, mov_file]);
mov_data = fread(mov_data_temp, 'single');
mov_frame_ind = mov_data([1:2:end]);
mov_data = mov_data([2:2:end]);

%alighn to cue onset
%plot cue aligned movement traces for each animal

%Align to first lick after the cue

%extract frame numbers corresponding to movement onset during trial

%align to movement onset during ITI

%compare movement amounts on imaging days to non-imaging days