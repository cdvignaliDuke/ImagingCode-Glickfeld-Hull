clear
file_info_CRP;
usFacs = 100;
behav_dir = 'Z:\Data\2P_imaging\behavior\';
for sub = [2, 5, 9, 22, 23, 31, 3, 6, 10, 25, 27, 32]%size(mouseID,2)
    for rID = 1
        file_info_CRP
        usFacs = 100; % upsample factor
        out_dir  = fullfile('Z:', 'Analysis','Cue_reward_pairing_analysis','2P',[dates{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
        
        data_dir = fullfile('Z:\Data\2P_imaging',[dates{sub} '_' mouseID{sub}], mouseID{sub});
        config_fn = dir(fullfile(data_dir,['*' runID{rID} '.mat']));
        img_fn   = dir(fullfile(data_dir,['*' runID{rID} '.sbx']));
        if size(img_fn,1) ~= 0
            [~,img_fn,~] = fileparts(img_fn.name);
            skip_run = 0;
            cd(data_dir);
            load(config_fn.name);
            nframes = info.config.frames;
        else
            skip_run = 1;
            data = [];
        end
        
        if skip_run == 1
            continue
        else
            if ~exist(out_dir)
                mkdir(out_dir);
            end
        end
        load([out_dir, '_cell_TCs']);
        load([out_dir, '_dFOverF_TC']);
        load([out_dir, 'ROI_TCs']);
        
        mouse = mouseID{sub};
        dateID = dates;
        dates = dateID{sub};
        
        subMat = dir([behav_dir '*' '9' mouse(end-1:end) '*' dates '*']);
        
        load([behav_dir, subMat.name]);
        dest = out_dir;
        
        getTC_events;
        CuePair_2P_TC_quantification;
        close all
        
    end
end