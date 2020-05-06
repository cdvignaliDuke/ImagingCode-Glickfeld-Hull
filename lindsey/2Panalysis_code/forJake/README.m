%% Below are the major analyses performed for the manuscript "Classical conditioning drives learned reward prediction signals in climbing fibers across the lateral cerebellum"
% Several functions used in these analysis are not listed, but are part of
% the ImagingCode-Glickfeld-Hull or BehaviorCode-Glickfeld-Hull
% repositories in the "Hull and Glickfeld Laboratories" organiation's
% GitHub page.

% % experiment lists:
% CRP_expt_list_all.m - datasets used for all main figures
% CRP_expt_list_Crus2019.m and CRP_expt_list_LS2019.m - datasets used for movement analysis

% % data processing (each dataset analyzed separately):
% CRP_extractTCs.m - motion registration, dendrite segmentation, timecourse extraction
% CRP_spikeAlign.m - conversion of dF/F into spikes and alignment to task events
% CRP_lickAlign.m - alignment of licking to task events
% CRP_lickResp.m - detection of spiking in response to lick events
% CRP_lickResp_preCueOnly - detection of spiking in response to lick events during ITI
% CRP_piezoAlign.m - alignment of movement signals to task events
% CRP_splitImage.m - GUI for segmenting lobules within FOV
% CRP_D1vsD2align.m - identification of same dendrites in pre- and post-learning sessions
% 
% % summary and figure making
% CRP_psthSummary.m - summary analysis of all spiking and licking data across all areas
% CRP_piezoSummary.m - summary analysis of spiking and movement data across all areas
% CRP_D1vsD2_Summary.m - comparison of spiking pre- vs post-learning for all neurons
% CRP_D1vsD2_Align_Summary.m - comparison of spiking pre- vs post-learning for same neurons
% CRP_example_cells.m - event aligned spiking for all trials for all cells
