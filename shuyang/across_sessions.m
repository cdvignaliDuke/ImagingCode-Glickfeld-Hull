%% SECTION ONE - assign pathnames and datasets to be analyzed/written. 
clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'180405_img1002_1'}; 
days = '1002-180405';
bx_source     = ['Z:\Data\Behv_MovingDots\behavior_raw'];
image_source_base  = ['Z:\Data\WF imaging\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF_MovingDots_Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the movingDots analysis folder
% behavior analysis results 
behav_output = ['Z:\Analysis\WF_MovingDots_Analysis\behavioral_analysis\' days];
color_code = {'b','r','k'};

%% SECTION TWO - draw df/f vs. speed across sessions for each individual


%% SECTION THREE - draw df/f vs. speed across all WTs and PCP2s


%% SECTION FIVE - df/f for each behavioral state across sessions for each individual


%% SECTION SIX - df/f for each behavioral state across all WTs and PCP2s