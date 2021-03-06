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
refData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [ref_date '_' mouse], [ref_date '_' mouse '_' ref_str], [ref_date '_' mouse '_' ref_str '_stimActFOV.mat']));
regData = load(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_RegData.mat']));
reg = regData.reg;
ref = regData.ref;
data_reg_to_ref = regData.data_reg_to_ref;
data_dfof_max_reg = regData.dfof_reg2ref;
data_dfof_max_ref = refData.data_dfof_max;
sz = size(reg);
rgb = zeros(sz(1),sz(2),3);

rgb(:,:,1) = data_dfof_max_ref;
rgb(:,:,2) = data_dfof_max_reg;
figure; image(rgb)
title([mouse ' ' date ' on ' ref_date ' RGB dfof overlay'])
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Reg2Ref_dfof.pdf']), '-dpdf')

rgb2 = zeros(sz(1),sz(2),3);
rgb2(:,:,1) = ref;
rgb2(:,:,2) = data_reg_to_ref;
figure; image(rgb)
title([mouse ' ' date ' on ' ref_date ' RGB data_reg overlay'])
print(fullfile('Z:\All_staff\home\grace\Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Reg2Ref_datareg.pdf']), '-dpdf')

end