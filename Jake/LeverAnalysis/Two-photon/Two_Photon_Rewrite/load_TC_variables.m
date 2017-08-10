%short script for loading the correct TC variables for use in
%peak_latency_2P and cell_resp_comparison

session_date = days_to_load{ii}(1:6);
session_mouseID = days_to_load{ii}(end-4:end);
session_runID = runID{1};
this_movie_dir = [TC_dir, session_date, '_', session_runID, '_', session_mouseID, '\_cue_movies.mat'];
this_category_dir = [TC_dir, session_date, '_', session_runID, '_', session_mouseID, '\_cell_categories.mat'];
this_TC_dir = [TC_dir, session_date, '_', session_runID, '_', session_mouseID, '\_cell_TCs.mat'];
this_lick_dir = [TC_dir, session_date, '_', session_runID, '_', session_mouseID, '\_cue_movies_lick.mat'];
if exist(this_movie_dir, 'file');
    load(this_movie_dir);
    load(this_category_dir);
    load(this_TC_dir);
    load(this_lick_dir);
else
    session_runID = runID{2};
    this_movie_dir = [TC_dir, session_date, '_', session_runID, '_', session_mouseID, '\_cue_movies.mat'];
    this_category_dir = [TC_dir, session_date, '_', session_runID, '_', session_mouseID, '\_cell_categories.mat'];
    this_TC_dir = [TC_dir, session_date, '_', session_runID, '_', session_mouseID, '\_cell_TCs.mat'];
    this_lick_dir = [TC_dir, session_date, '_', session_runID, '_', session_mouseID, '\_cue_movies_lick.mat'];
    load(this_movie_dir);
    load(this_category_dir);
    load(this_TC_dir);
    load(this_lick_dir);
end