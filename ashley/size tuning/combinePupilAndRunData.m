clear all
close all
ds = 'szTuning_PM';
eval(ds)

rc = behavConstsAV;
% imgParams_FSAV
fn = 'Z:\home\ashley\Manuscripts\Size Tuning\Matlab Analysis';

load(fullfile(fn,'eyeSummary'))
load(fullfile(fn,'runSummary'))

nexp = size(expt,2);
eyeRun = struct;
for iexp = 1:nexp
    eyeRun(iexp).mouse = eyeExpt(iexp).mouse;
    eyeRun(iexp).date = eyeExpt(iexp).date;
    eyeRun(iexp).sz.tSz = eyeExpt(iexp).tSz;
    eyeRun(iexp).sz.relEyePosition = eyeExpt(iexp).pos(3).position;
    eyeRun(iexp).sz.speedCPS = runExpt_cells(iexp).sz.allTrialsSpeedCPS;
    eyeRun(iexp).ret.tAz = eyeExpt_ret(iexp).tAz;
    eyeRun(iexp).ret.tEl = eyeExpt_ret(iexp).tEl;
    eyeRun(iexp).ret.relEyePosition = eyeExpt_ret(iexp).pos(3).position;
    eyeRun(iexp).ret.speedCPS = runExpt_cells(iexp).ret.allTrialsSpeedCPS;
end

save(fullfile(fn,'eye&RunStruct'),'eyeRun')