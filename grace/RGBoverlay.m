%% clear everything
clear all
clear global
%% get path names
date = '200108';
ImgFolder = strvcat('003');
run = strvcat('000');
mouse = 'i1316';
doFromRef = 1;
ref_date = '200106';
ref_run = strvcat('003');
nrun = size(ImgFolder,1);
frame_rate = 15.5;
run_str = catRunName(ImgFolder, nrun);
%% RGB Overlay
if doFromRef
    ref_str = ['runs-' ref_run];
    if size(ref_run,1)>1
        ref_str = [ref_str '-' ref_run(size(ref_run,1),:)];
    end
stimActFOVref = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']));
stimActFOVreg = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_dfof_max.mat']));
data_dfof_max_ref = stimActFOVref.data_dfof_max;
data_dfof_max_reg = stimActFOVreg.dfof_reg2ref;
data_dfof_avg_all = stimActFOVref.data_dfof_avg_all;
sz = size(data_dfof_avg_all);
rgb = zeros(sz(1),sz(2),3);
rgb(:,:,1) = data_dfof_max_ref;
rgb(:,:,2) = data_dfof_max_reg;
figure; image(rgb)
title([mouse ' ' date ' on ' ref_date ' RGB overlay'])
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_RegOnRefFOV.pdf']), '-dpdf')
end