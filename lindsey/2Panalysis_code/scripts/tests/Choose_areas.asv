clear all

date = '111125';
mouse = 'M14';
userun = [3:5];

base = 'G:\users\lindsey\analysisLG\active mice';
running_base = 'G:\users\lindsey\dataLG\Running data';
outDir = fullfile(base, mouse,date);
%% load image for choosing rois
fn = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_max_blur.tif']);
max_dF = readtiff(fn);
bwimgarea = imCellEditInteractive2(max_dF,[]);
area_mask = bwlabel(bwimgarea);
figure; imagesc(area_mask);

a = 'LM'; b = 'AL'; c = 'V1'; d = 'RL';  e = 'A'; f = 'AM'; g = 'PM';
area_list = strvcat(a,b,c,d,e,f,g);

fn_out = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_area_mask.mat']);
save(fn_out, 'area_mask', 'area_list');
