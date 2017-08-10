clear all
close all
%%
doDendrites = 0;
rc = behavConstsAV;
awFSAVdatasets_V1gad
slct_expt = 1:size(expt,2);
%%
for iexp = slct_expt
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
if any(strcmp(fieldnames(expt),'nTunFrames'))
    [input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName,expt(iexp).nTunFrames);  
else
    [input, data] = Load_SBXdataPlusMWorksData(SubNum,expDate,dirTime,mouse,dirFolder,fName);
end
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
if doDendrites
    load(fullfile(fn,'dendrite_mask.mat'))
else
    load(fullfile(fn,'final_mask.mat'))
end
%% get timecourses and subtract neuropil
buf = 4;
np = 6;
if strcmp(expt(iexp).img_strct,'axons')
    data_tc = stackGetTimeCourses(data_reg,mask_boutons);  
    data_tc_subnp = getWeightedNeuropilTimeCourse(data_reg,data_tc,mask_boutons,buf,np);  
else
    data_tc = stackGetTimeCourses(data_reg,mask_cell);
    data_tc_subnp = getWeightedNeuropilTimeCourse(data_reg,data_tc,mask_cell,buf,np);
end

%% save data
if doDendrites
    save(fullfile(fntun,'timecourses_dendrites.mat'),'data_tc_subnp')
    save(fullfile(fntun,'raw_tc_dendrites.mat'),'out_tun','data_tc','buf','np')
else
    save(fullfile(fntun,'timecourses.mat'),'data_tc_subnp')
    save(fullfile(fntun,'raw_tc.mat'),'out_tun','data_tc','buf','np')
end
quickDirectionTuningCurves
end