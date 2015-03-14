%% animal and experiment date
SubNum = '516';
mouse = '516';
date = '141003';

%% open notes pdf
CD = ['Z:\data\' mouse '\two-photon imaging\' date];
cd(CD);
notes = ['2P Notes - ' SubNum ' ' date '.pdf'];
open(notes);

%% Retinotopy

%load retinotopy dataset
edit Load_SBXdataset_fast.m

%find preferred retinotopy - run up to section 'response matrix for each
%position'
edit RetinotopyAnalysis.m

%save txt notes for mat files

%% Direction Tuning

%load direction tuning dataset
edit Load_SBXdataset_fast.m

%create mask that will be used for finding cells in flashing stim exp
%save averaged image to register flashing stim datasets to
%find tuning properties of those cells

edit DirectionTuning_V1.m

%plot response to each direction for all cells

edit DirectionTuning_V1_plotAllDirectionsAllCells.m

%plot population response to each direction

edit DirectionTuning_V1_population.m

%% Flashing Stim Analysis
