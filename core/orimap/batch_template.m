% Caution! Create a new directory & set the current directory there

Input_fname='C:\Share\OG031210\OG031210_3\OG031210_3_vstim8_t';

nframes_per_stim=10;        % number of frames per one stimulation
nstim_per_run=8;            % 8 directions of gratings
nframes_per_trial = nframes_per_stim .* nstim_per_run;

run_inds=[0:3,5:9];         % run indexes to be averaged
nbinning=2;                 % spatial binning factor

base_inds=[1:4];            % baseline frames in one trial (between 1 and nframes_per_stim)
stim_inds=[5:9];            % activation frames in one trial (between 1 and nframes_per_stim)

lowpass = 3;
kernel_size = 15;

sp_filter=fspecial('gaussian', kernel_size, lowpass);    % spatial lowpass filter in pixel unit (after binning)

tune_max = 0.4; % in HLS map, if tune >= tune_max, color will be saturated. if not, color will be unsaturated. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ave=read_stack_ave(Input_fname,nframes_per_trial,run_inds,nbinning);
save ave ave;

[base, base_sm] = calc_base (ave, base_inds, nframes_per_stim, nstim_per_run, sp_filter);
save base base;
save base_sm base_sm;

dir_F=calc_F(ave, stim_inds, nframes_per_stim, nstim_per_run);
save dir_F dir_F;

[dir_dF, dir_dF_sm, ori_dF, ori_dF_sm] = calc_dF (dir_F, base, sp_filter);
save dir_dF dir_dF;
save dir_dF_sm dir_dF_sm;
save ori_dF ori_dF;
save ori_dF_sm ori_dF_sm;

[dir_ratio, ori_ratio] = calc_ratio_change (dir_dF_sm, ori_dF_sm, base_sm);
save dir_ratio dir_ratio;
save ori_ratio ori_ratio;

[max_dir_dF, max_dir_ratio_change, max_ori_dF, max_ori_ratio_change] = write_dF_images (dir_dF_sm, ori_dF_sm, dir_ratio, ori_ratio);

dir_dF_params = calc_map_params (dir_dF_sm);
dir_ratio_params = calc_map_params (dir_ratio);

ori_dF_params = calc_map_params (ori_dF_sm);
ori_ratio_params = calc_map_params (ori_ratio);

save dir_dF_params dir_dF_params;
save dir_ratio_params dir_ratio_params;
save ori_dF_params ori_dF_params;
save ori_ratio_params ori_ratio_params;

write_intensity_maps (dir_dF_params, 'dir_dF_', max_dir_dF);
write_intensity_maps (dir_ratio_params, 'dir_ratio_', max_dir_ratio_change);
write_intensity_maps (ori_dF_params, 'ori_dF_', max_ori_dF);
write_intensity_maps (ori_ratio_params, 'ori_ratio_', max_ori_ratio_change);

dir_angle = write_angle_map (dir_dF_params.th, 'dir_');
ori_angle = write_angle_map (ori_dF_params.th, 'ori_');

save dir_angle dir_angle;
save ori_angle ori_angle;

dir_dF_polar = write_polar_map (dir_dF_params, 'dir_dF_', max_dir_dF);
dir_ratio_polar = write_polar_map (dir_ratio_params, 'dir_ratio_', max_dir_ratio_change);
ori_dF_polar = write_polar_map (ori_dF_params, 'ori_dF_', max_ori_dF);
ori_ratio_polar = write_polar_map (ori_ratio_params, 'ori_ratio_', max_ori_ratio_change);

save dir_dF_polar dir_dF_polar;
save dir_ratio_polar dir_ratio_polar;
save ori_dF_polar ori_dF_polar;
save ori_ratio_polar ori_ratio_polar;

dir_dF_HLS = write_HLS_map (dir_dF_params, 'dir_dF_', max_dir_dF, tune_max);
dir_ratio_HLS = write_HLS_map (dir_ratio_params, 'dir_ratio_', max_dir_ratio_change, tune_max);
ori_dF_HLS = write_HLS_map (ori_dF_params, 'ori_dF_', max_ori_dF, tune_max);
ori_ratio_HLS = write_HLS_map (ori_ratio_params, 'ori_ratio_', max_ori_ratio_change, tune_max);

save dir_dF_HLS dir_dF_HLS;
save dir_ratio_HLS dir_ratio_HLS;
save ori_dF_HLS ori_dF_HLS;
save ori_ratio_HLS ori_ratio_HLS;

log.input_file = Input_fname;
log.nframes_per_stim = nframes_per_stim;
log.nstim_per_run = nstim_per_run;
log.nframes_per_trial = nframes_per_trial;
log.run_inds = run_inds;
log.nbinning = nbinning;
log.base_inds = base_inds;
log.stim_inds = stim_inds;
log.lowpass = lowpass;
log.kernel_size = kernel_size;
log.tune_max = tune_max;
log.max_dir_dF = max_dir_dF;
log.max_dir_ratio_change = max_dir_ratio_change;
log.max_ori_dF = max_ori_dF;
log.max_ori_ratio_change = max_ori_ratio_change;

save log log;
