clear all
close all
ds = 'awFSAVdatasets_V1gad';
rc = behavConstsAV;
eval(ds)
iexp = 8
% for iexp = slct_exp
%%
mouse = expt(iexp).mouse;
expDate = expt(iexp).date;
redFolder = expt(iexp).redChannelRun;
redChannelOn = expt(iexp).greenredsimultaneous;
disp([mouse expDate redFolder])
fn = fullfile(rc.ashleyAnalysis,mouse,'two-photon imaging', expDate,'data processing');

%% load red data
fName = [redFolder '_000_000'];
redData = loadsbx_choosepmt(2,mouse,expDate,redFolder,fName,1000);
figure;imagesc(mean(redData,3));colormap gray
% end
greenData = loadsbx_choosepmt(1,mouse,expDate,redFolder,fName,1000);

%% register
load(fullfile(fn,'regOuts&Img.mat'))
if redChannelOn == 1
    [~,redReg] = stackRegister(redData,regImg);
else
    [outs,greenReg] = stackRegister(greenData,regImg);
    outs = double(outs);
    [~,redReg] = stackRegister_MA(redData,[],[],outs);
end
clear greenData redData
redImage = mean(redReg,3);
figure;colormap hot; imagesc(redImage)
figure;colormap hot; imagesc(regImg)
%% save image
save(fullfile(fn,'redChannelImage'),'redImage')
% end