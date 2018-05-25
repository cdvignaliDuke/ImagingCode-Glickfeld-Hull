clear all
close all
ds = 'FSAV_V1_GAD';
rc = behavConstsAV;
eval(ds)
for iexp = 1:size(expt,2)
%%
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
redFolder = expt(iexp).redChannelRun;
redChannelOn = expt(iexp).greenredsimultaneous;
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,'data processing');

%% load red data
fName = [redFolder '_000_000'];
redData = loadsbx_choosepmt(2,mouse,expDate,redFolder,fName,1000);

greenData = loadsbx_choosepmt(1,mouse,expDate,redFolder,fName,1000);

%% register
load(fullfile(fn,'regOuts&Img.mat'))
if redChannelOn == 1
    [~,redReg] = stackRegister(redData,regImg);
else
    [outs,greenReg] = stackRegister(greenData,regImg);
    [~,redReg] = stackRegister_MA(redData,[],[],outs);
end
clear greenData redData
redImage = mean(redReg,3);
figure;colormap hot; imagesc(redImage)

%% save image
save(fullfile(fn,'redChannelImage'),'redImage')
end