%% write frames around running start into a movie,
% to test if the image shifts when the mice starts moving
% UNCOMPLETED!

%% file paths
clear;
sessions = '190422_img1015'; 
image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\'; %stores the data on crash in the movingDots analysis folder
image_analysis_dest = [image_analysis_base, sessions, '\'];

% behavior analysis results 
days = '1015-190422_1';
behav_dest = ['Z:\Analysis\2P_MovingDots_Analysis\behavioral_analysis\' days '\'];
color_code = {'b','r','k','c'};


%%
behav_output = load([behav_dest days '_behavAnalysis.mat']);
frm_run_cell = behav_output.frames_run_cell;
frm_run = cell2mat(frm_run_cell);
runonset_inx = find(diff(frm_run)>1)+1;
runonest = frm_run(runonset_inx);

cd(image_source);
order = '000';
%file = [ID '_000_' order];
file = [sessions '_000_' order];
frame_start = 0;
nframes = 29905;
imgread = squeeze(sbxread(file,frame_start,nframes));
f1 = 1;
f2 = 1000;
img_tiff = imgread(:,:,f1:f2);
writetiff(img_tiff,[image_analysis_dest sessions '_'...
    order '_tiff_', num2str(f1), '_', num2str(f2)]);


writetiff(img_rgs(:,:,idx),[image_analysis_dest sessions '_'...
    order '_rgstr_tiff_', num2str(frame_start), '_', num2str(nframes) '_',num2str(gap) '_ref' num2str(ref_frame)]);
save([image_analysis_dest sessions '_' order '_img_rgs_ref' num2str(ref_frame) '.mat'], 'img_rgs', 'rgs_out', 'img_ref','-append');




%% Nathan's code
% Find relevant behavioral event to align to and extract frames around it
mworks = input; clear input
sz = size(img_reg);
cTargetOn = celleqel2mat_padded(mworks.cTargetOn);
nTrials = length(cTargetOn);
prewin_frames = round(1500/mworks.frameRateHz);
postwin_frames = round(3000/mworks.frameRateHz);
reg_align = nan(sz(1),sz(2),prewin_frames+postwin_frames,nTrials);
for itrial = 1:nTrials
    if cTargetOn(itrial)+postwin_frames < size(img_reg,3)
        reg_align(:,:,:,itrial) =img_reg(:,:,cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1);
    end
end
reg_align_avg = nanmean(reg_align,4);
writetiff(reg_align_avg,fullfile(outdir,img_fn,['Registered-Movie-',register_type,'.tif']));
fprintf('done \n');