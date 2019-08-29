clear all
CRP_expt_list
nd = size(expt,2);
behav_dir = '\\crash.dhe.duke.edu\data\home\andrew\Behavior\Data';
jake_dir = '\\crash.dhe.duke.edu\data\home\jake\Analysis\Cue_reward_pairing_analysis\2P';
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';

for id =4
    fprintf(['Day: ' num2str(id) '\n'])
    nexp = size(expt(id).date,1);
    [b b2] = subplotn(nexp);
    figure;
    for iexp = 1:nexp
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([mouse ' ' date '\n'])
        img_fn = [date '_' mouse];
        
        if str2num(mouse)>800
            mouse_name = mouse(2:3);
        else
            mouse_name = mouse;
        end
        
%         cd(fullfile('\\crash.dhe.duke.edu\data\home\jake\Data\2P_imaging', [date '_img' mouse_name],['img' mouse_name]))
%         load(['img' mouse_name '_000_' run '.mat'])
%         load(fullfile(lg_out, img_fn, [img_fn '_reg.mat']));
%         load(fullfile(lg_out, img_fn, [img_fn '_input.mat'])); 
%         sz = size(img_ref);    
%         cTargetOn = celleqel2mat_padded(input.cTargetOn);
%         nTrials = length(cTargetOn);
%         prewin_frames = round(1500./input.frameRateHz);
%         postwin_frames = round(3000./input.frameRateHz);
%         reg_align = nan(sz(1),sz(2),prewin_frames+postwin_frames,nTrials);
%         for itrial = 1:nTrials
%             if cTargetOn(itrial)+postwin_frames<input.counterValues{end}(end) & cTargetOn(itrial)+postwin_frames<input.counterValues{itrial}(end)
%                 data = sbxread(['img' mouse_name '_000_' run],cTargetOn(itrial)-prewin_frames,postwin_frames+prewin_frames);
%                 data = squeeze(data);
%                 [out img_reg] = stackRegister_MA(data,[],[],outs(cTargetOn(itrial)-prewin_frames:cTargetOn(itrial)+postwin_frames-1,:));
%                 reg_align(:,:,:,itrial) = img_reg;
%             end
%         end
%         reg_align_avg = nanmean(reg_align,4);
%         writetiff(reg_align_avg(:,:,1:2:end), fullfile('\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_movies', ['Day' num2str(id)], [img_fn '_cueAlign_rereg.tif']))
%         

        load(fullfile(lg_out, img_fn, [img_fn '_ROI_TCs.mat']));
        figure; imagesc(mask_flat); truesize
        title([mouse ' ' date])
        print(fullfile('\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_mask', ['Day' num2str(id)], [img_fn '_finalMask.pdf']),'-dpdf', '-fillpage')

        openfig(fullfile(lg_out,img_fn, [img_fn '_cueAlign_events_Hz.fig']));
        print(fullfile('\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_spikeDetect', ['Day' num2str(id)], [img_fn '_cueAlign_events_Hz.pdf']),'-dpdf', '-fillpage')

        openfig(fullfile(lg_out,img_fn, [img_fn '_cueAlign_lickHz.fig']));
        print(fullfile('\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_licking', ['Day' num2str(id)], [img_fn '_cueAlign_lickHz.pdf']),'-dpdf', '-fillpage')
        
        close all
    end
end
    
