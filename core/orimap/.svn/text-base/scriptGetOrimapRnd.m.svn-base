%exptName='mouse090702';
%seriesName='vstim3';

load rndSeq090522

%%%%%%%%%%%%%%%%%%%
% stim parameters
%%%%%%%%%%%%%%%%%%%

Decimation = 30;   % temporal decimation

Noff = 4;
Non = 4;
nstim_per_run=12;            % 8 directions of gratings

discard_off_frames = 2; % How many off frames to discard to avoid calcium tail
                        % recommendtation 2sec (convert it to frame number)

%%%%%%%%%%%%%%%%%%%%                        
% map parameters
%%%%%%%%%%%%%%%%%%%
                        
nbinning=1;                 % spatial binning factor

lowpass = 0.5; %default=1.5
kernel_size = 3;
sp_filter=fspecial('gaussian', kernel_size, lowpass);    % spatial lowpass filter in pixel unit (after binning)

tune_max = 0.4; % in HLS map, if tune >= tune_max, color will be saturated. if not, color will be unsaturated. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% preprocessing of data for fastrig
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% decimate, if not on disk already
dirs = frGetDirs;
decDir = fullfileforw(dirs.images, exptName, seriesName, ...
                      ['binned_green']);
decFile = fullfileforw(decDir, ['binned_' seriesName '.tif']);
if ~exist(decFile, 'file')
    [stack, outdir, fname] = frReadWithDecimation(exptName,  seriesName, ...
                                                  'green', true, Decimation);
else
    stack = readtiff(decFile);
end
target=sum(stack,3);

% register, if not on disk already
regFile = fullfileforw(decDir, ['rbinned_' seriesName '.tif']);
if ~exist(regFile, 'file')
    [realignParams,rstack]=stackRegister(stack,target);
    writetiff(rstack, regFile);
    save(fullfileforw(decDir, 'realignParams'), 'realignParams');
else
    rstack = readtiff(regFile);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reverse correlation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


delays=[discard_off_frames:Noff*2+Non-1];
                % from the middle of prestim off period to the end of
                % poststim off period
rcimg=stackRevCor(rstack,condition,delays,Noff+Non);
writetiff(rcimg,fullfileforw(decDir,'rcimg.tif'));

nframes_per_stim=length(delays);        % number of frames per one stimulation
nframes_per_trial = nframes_per_stim .* nstim_per_run;


run_inds=[0];         % run indexes to be averaged

base_inds=[1:Noff-discard_off_frames];            % baseline frames in one trial (between 1 and nframes_per_stim)
stim_inds=[Noff-discard_off_frames+1:Noff-discard_off_frames+Non]; % activation frames in one trial (between 1 and nframes_per_stim)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

cd(outdir);

mkdir(outdir, 'maps');
cd maps

[Nx,Ny,Nt]=size(rcimg)
ave=double(rcimg);
save ave ave;

FOV = sum(ave,3);
save FOV FOV;
imwrite(FOV./max(max(FOV)), 'FOV.tif');

[Nx,Ny,Nt]=size(ave);
bases=calc_F(ave, base_inds, nframes_per_stim, nstim_per_run);
bases_sm=zeros(Nx,Ny,nstim_per_run);
for i=1:nstim_per_run
    bases_sm(:,:,i)=imFilter2(sp_filter,bases(:,:,i));
end

save bases bases;
save bases_sm bases_sm;

dir_F=calc_F(ave, stim_inds, nframes_per_stim, nstim_per_run);
for i=1:nstim_per_run
    bases_sm(:,:,i)=imFilter2(sp_filter,bases(:,:,i));
end

save dir_F dir_F;


[dir_dF, dir_dF_sm, ori_dF, ori_dF_sm] = calc_dFb (dir_F, bases, sp_filter);
save dir_dF dir_dF;
save dir_dF_sm dir_dF_sm;
save ori_dF ori_dF;
save ori_dF_sm ori_dF_sm;

[dir_ratio, ori_ratio] = calc_ratio_changeb (dir_dF_sm, ori_dF_sm, bases_sm);
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

log.fname=fname;
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
