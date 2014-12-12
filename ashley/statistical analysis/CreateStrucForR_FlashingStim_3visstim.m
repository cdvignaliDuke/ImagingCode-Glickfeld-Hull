%% structure variables
% cell == the structure name and cell(1) indicated data relevant to that
% cell, possibly accross multiple experiments, but obviously only one mouse
% animal == animal ID number
% training == naive mouse (naive) or mouse that has learned the task
% (trained)
% auditory == auditory context, vector containing dF/F value for 42 frames of visual
% stimulus presentation 
% visual == visual context, vector containing dF/F value for 42 frames of visual
% stimulus presentation
% oriPref == preferred orientation of that cell, an integer
% oriOSI == orientation selectivity index of that cell
%% 
%load bx analysis mat file
CD = 'Z:\2P imaging\Analysis\516\141003\002+003+004\FlashingStimAnalysis\Rsp2VisStim';
cd(CD);
load('analysis.mat');

% *build training portion of structure*
%load dF/F vectors
for icell = 1:nCells
    cell(icell).auditory = Cell_3visstim_mean_A(:,icell);
end

for icell = 1:nCells
    cell(icell).visual = Cell_3visstim_mean_V(:,icell);
end
% *build orientation selectivity of the structure*
% load direction selectivity analysis
CD = 'Z:\2P imaging\Analysis\516\141003\005 - Direction Selectivity'
cd(CD);
load('analysis.mat');

for icell = 1:nCells
    cell(icell).oriPref = CellPref_deg(icell);
end

for icell = 1:nCells
    cell(icell).oriOSI = CellOSI(icell);
end


for icell = 1:nCells
    cell(icell).animal = '516';
end

for icell = 1:nCells
    cell(icell).training = 'trained';
end

CD = 'Z:\2P imaging\Analysis\Data Structures\141013';
cd(CD);
save('cellStruct.mat','cell');

%% **control mouse**

%load flashing stim analysis
nCells = size(cell,2);
tCells = size(data_TC,2);
nCells = nCells+1:nCells+tCells;

start = 1;
for icell = nCells
    cell(icell).auditory = Cell_3visstim_mean_A(:,start);
    start = start+1;
end

start = 1;
for icell = nCells
    cell(icell).visual = Cell_3visstim_mean_V(:,start);
    start = start+1;
end

%load direction selelctivity analysis
CD = 'Z:\2P imaging\Analysis\AW04\140923\003 - Direction Selectivity';
cd(CD);
load('analysis.mat');


start = 1;
for icell = 1:nCells
    cell(icell).oriPref = CellPref_deg(icell);
    start = 1;
end

start = 1;
for icell = 1:nCells
    cell(icell).oriOSI = CellOSI(icell);
    start = 1;
end

%%
%****make excel sheet****
nCells = size(cell,2);
Cell = 1:nCells;
animal = {cell.animal};
training = {cell.training};
oriPref = {cell.oriPref};
OSI = {cell.oriOSI};
for icell = Cell
    auditory(:,icell) = cell(icell).auditory;
end
for icell = Cell
    visual(:,icell) = cell(icell).visual;
end

%aud 1-10
for icell = Cell
    A1(icell) = mean(auditory(1:10,icell),1);
end
%11-21
for icell = Cell
    A2(icell) = mean(auditory(11:21,icell),1);
end
%22-31
for icell = Cell
    A3(icell) = mean(auditory(22:31,icell),1);
end
%32-42
for icell = Cell
    A4(icell) = mean(auditory(32:42,icell),1);
end
%vis 1-10
for icell = Cell
    V1(icell) = mean(visual(1:10,icell),1);
end
%11-21
for icell = Cell
    V2(icell) = mean(visual(11:21,icell),1);
end
%22-31
for icell = Cell
    V3(icell) = mean(visual(22:31,icell),1);
end
%32-42
for icell = Cell
    V4(icell) = mean(visual(32:42,icell),1);
end

for icell = Cell
    cell(icell).A1 = A1(icell);
    cell(icell).A2 = A2(icell);
    cell(icell).A3 = A3(icell);
    cell(icell).A4 = A4(icell);
    cell(icell).V1 = V1(icell);
    cell(icell).V2 = V2(icell);
    cell(icell).V3 = V3(icell);
    cell(icell).V4 = V4(icell);
end

x = struct2table(cell);

filename = 'cellStruct.xlsx';
writetable(x,filename,'Sheet',1,'Range','D1')