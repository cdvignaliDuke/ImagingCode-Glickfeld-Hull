clear all
close all
rc = behavConstsAV;
awFSAVdatasets_V1
for iexp = [7,10,12,13,18,19,20]
% iexp = 11;

SubNum = expt(iexp).SubNum;
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
dirFolder = expt(iexp).dirtuning;
dirTime = expt(iexp).dirtuning_time;
down = 10;

fntun = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate, dirFolder);
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate);
%% load tuning data
fName = [dirFolder '_000_000'];
[input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName);  

% down-sample
data_down = stackGroupProject(data,down);
clear data

% remove negative data by subtraction
data_sub = data_down-min(min(min(data_down,[],1),[],2),[],3);
clear data_down


%% load outs and re-regester data
load(fullfile(fn,'regOuts&Img.mat'));

[out data_reg] = stackRegister_MA(data_sub,[],[],out_tun);

%% load mask
load(fullfile(fn,'final_mask.mat'))
%% get timecourses
data_tc = stackGetTimeCourses(data_reg,mask_cell);

%% subtract neuropil

% get neuropil timecourses
buf = 4;
np = 6;

data_tc_subnp = getWeightedNeuropilTimeCourse(data_reg,data_tc,mask_cell,buf,np);

%% save data
save(fullfile(fntun,'timecourses.mat'),'data_tc_subnp')
save(fullfile(fntun,'raw_tc.mat'),'out_tun','data_tc','buf','np')

quickDirectionTuningCurves
end