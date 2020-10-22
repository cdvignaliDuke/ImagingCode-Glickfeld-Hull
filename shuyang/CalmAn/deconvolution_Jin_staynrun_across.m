% spike identification using deconvolution. (function deconvolve Ca got from OASIS on git hub)
% to use the deconvolution function, need to add OASIS function to search
% path: >>oasis_setup
% function deconvolve Ca: denoise raw fluorescence trace and identifies spike events. 
% right now using FOOPSI method and the model is ar1 (also tried ar2: see commented bottom part),ar1 works better than ar2.
% then SECTION II draws the GUI showing raw traces with identified events and deconvolved trace.
% 
%% assign document paths and experimental sessions

clear;
sessions = {'190429_img1021','190430_img1023','190507_img1024','190603_img1025'};
days = {'1021-190429_1','1023-190430_1','1024-190507_1','1025-190603_1'};

image_analysis_base    = 'Z:\Analysis\2P_MovingDots_Analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder


%%
btm_peak_amp_stay_across = zeros(1,length(sessions));
btm_peak_amp_run_across = zeros(1,length(sessions));
for ii = 1:length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    image_analysis_dest_deconv_sep = [image_analysis_base, sessions{ii}, '\deconvolution_sep\'];
    decon_sep_output = load([image_analysis_dest_deconv_sep sessions{ii} '_spk_deconvolve_staynrun_seperate.mat' ]);
    peak_amp_stay = decon_sep_output.peak_amp_stay;
    peak_amp_run = decon_sep_output.peak_amp_run;
    btm_peak_amp_stay_across(ii) = peak_amp_stay(10);
    btm_peak_amp_run_across(ii) = peak_amp_run(10);
end

