clear all
close all
for id = [1:5]; 
   % CRP_OT_expt_list_LS_jh;
    CRP_OT_expt_list_Crus_jh;
    lg_out = 'Y:\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';
    behav_dir = 'Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior';
    img_dir = 'Y:\home\jake\Data\2P_imaging';
    nexp = size(expt(id).date,1);
    fprintf(['Day: ' num2str(id) '\n'])
    for iexp = 1:nexp
        %define session identifiers
        mouse = strtrim(expt(id).mouse(iexp,:));
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];
        
        %check to see if piezo data exists for this session
        if ~exist([img_dir, '\', [date '_img' mouse], '\', ['img' mouse], ['\img' mouse '_000_' run], '.ephys'])
            fprintf([img_fn, ' has no piezo data. Continuing to next session. \n' ]);
            continue
        end

        %set directory and load data
        cd(fullfile(img_dir, [date '_img' mouse],['img' mouse]))
        load(fullfile(lg_out,img_fn, [img_fn '_piezo_frames.mat']), 'piezo_frames');
        piezo_std = std(piezo_frames);
        piezo_mean = mean(piezo_frames);
        piezo_outlie_ind = find(piezo_frames > piezo_mean+(piezo_std*15) );
        
        %highpass filter
        piezo_hp_filt = highpass(piezo_frames, 1, 30);
        
%         %remove outliers
%         if ~isempty(piezo_outlie_ind)
%             piezo_frames(piezo_outlie_ind) = NaN;
%         end
%         piezo_z = zscore(piezo_frames);
         
        
        %% plot outputs for each session
        figure;
        subplot(1,2,1); 
        plot(piezo_frames); %hold on;
        title('Volts')
        ylabel('Volts');
%         
%         subplot(1,2,2);
%         plot(piezo_z);
%         title('Z score');
%         ylabel('Z score');
        
        subplot(1,2,2);
        plot(piezo_hp_filt);
        title('highpass filter: 1 Hz');
        ylabel('Volts');
        hline(piezo_std, 'r')
        suptitle([date ' ' mouse, ' LS piezo data']);
        
        
        
    end
end