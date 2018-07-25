% pointRegistration_meta script
%script for using the Thresholding method to extract ROIs from 2P CRP data

%% Set directories and load the data
clear
%file_info_CRP;
% use 180425_img085 180502_img085       180507_img081  180518_img081
dates = {'180425'};
mouseID = {'img084'};
behavID = {'084'};
runID = {'000', '001', '002'};
irun = 1;

behav_dir = 'Z:\Data\2P_imaging\behavior_CRP\';
crp_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\FoV_matching_outputs_thresholding';
docombineMask = 0; doRerun = 0; useWarp = 0; doPCA = 0;
useGPU = 0; checkImg = 0;
rID = 1; sub =1;

%% Motion Registration

%set pathnames
usFacs = 100; % upsample factor for stackRegister
data_dir = fullfile('Z:\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub}, '\');
out_dir =  fullfile(crp_dir, [dates{sub}, '_', runID{rID}, '_', mouseID{sub}], '\');
if ~exist(out_dir)
    mkdir(out_dir);
end

% Extract pockel cell data
if exist([data_dir, mouseID{sub}, '_', runID{rID}, '_', runID{rID}, '_realtime.mat'], 'file') ~= 2 %checks for data about pockel cell
    %set directories
    mouse = mouseID{sub};
    subMat = dir([behav_dir, '*', behavID{sub}, '-', dates{sub}, '*']); %sets directory as the behavior file
    
    %load the behavior data
    load([behav_dir, subMat.name]);
    
    %identify times of "lever press" and trial end
    cStart = cell2mat(input.cLeverDown);
    cEnd = cell2mat(input.cTrialEnd);
    
    %load imaging matfile
    config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
    load([data_dir, config_fn.name]);
    
    %for some reason Ziye overwrote ttl_log with this new version derived from the bx data. It makes no sense to I renamed it to ttl_log2 so it wont overwrite ttl_log
    ttl_log2 = zeros(info.config.frames,1);
    for jj = 1:length(cStart)
        ttl_log2(cStart(jj):cEnd(jj)) = 1;
    end
    %save([data_dir, mouse, '_', runID{rID}, '_realtime.mat'], 'ttl_log2');  %overwrites existing ttl_log data
end

%Motion registration or check for existant motion reg outputs
if exist([out_dir, 'Reg_out.mat'],'file') == 2
    load([out_dir, 'Reg_out.mat']);
    load([out_dir, 'img_reg.mat']);
    %[~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,double(reg_out));
elseif exist([out_dir, 'img_reg.mat'],'file') == 2
    load([out_dir, 'img_reg.mat']);
else
    [img_reg, laser_on_ind_conserv, nt] = motionCorrect(data_dir, runID{rID}, out_dir); %uses stackRegister, automatically saves img_reg and Reg_out
end


%% Threshold method of dendrite identification

%pointRegistration



