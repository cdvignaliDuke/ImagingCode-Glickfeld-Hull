clear all
close all
CRP_expt_list_Crus2019
for id = 2 
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
jake_dir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\jake\Analysis\Cue_reward_pairing_analysis\CC_analysis_2P';
nexp = size(expt(id).date,1);
fprintf(['Day ' num2str(id) '\n'])
    for iexp = 1:nexp
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];

        if exist(fullfile(lg_out,img_fn, [img_fn '_reg.mat']))
            load(fullfile(lg_out,img_fn, [img_fn '_reg.mat']))
        elseif exist(fullfile(jake_dir,img_fn, [img_fn '_reg.mat']))
            load(fullfile(jake_dir,img_fn, [img_fn '_reg.mat']))
        end
            
        sz = size(img_ref);
%         imagesc(img_ref)
        load(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))
        if size(expt(id).areas{iexp},2)>1
            if exist(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
                load(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            elseif exist(fullfile(jake_dir,img_fn, [img_fn '_splitImage.mat']))
                load(fullfile(jake_dir,img_fn, [img_fn '_splitImage.mat']))
            end
            figure;
            indL = find(maskCat==1);
            indR = find(maskCat==2);
            subplot(2,2,1)
            shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,indL,ind_rew),3),2), (nanstd(nanmean(targetAligndFoverF(:,indL,ind_rew),3),[],2))./sqrt(length(indL)),'k');
            xlabel('Time from cue')
            ylabel('dF/F')
            title(['Reward- Left side- n=' num2str(length(indL))])
            vline(600)
            subplot(2,2,2)
            shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,indR,ind_rew),3),2), (nanstd(nanmean(targetAligndFoverF(:,indR,ind_rew),3),[],2))./sqrt(length(indR)),'k');
            xlabel('Time from cue')
            ylabel('dF/F')
            vline(600)
            title(['Reward- Right side- n=' num2str(length(indR))])
            subplot(2,2,3)
            shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,indL,ind_omit),3),2), (nanstd(nanmean(targetAligndFoverF(:,indL,ind_omit),3),[],2))./sqrt(length(indL)),'r');
            xlabel('Time from cue')
            ylabel('dF/F')
            title(['Omit- Left side- n=' num2str(length(indL))])
            vline(600)
            subplot(2,2,4)
            shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,indR,ind_omit),3),2), (nanstd(nanmean(targetAligndFoverF(:,indR,ind_omit),3),[],2))./sqrt(length(indR)),'r');
            xlabel('Time from cue')
            ylabel('dF/F')
            title(['omit- Right side- n=' num2str(length(indR))])
            vline(600)
            suptitle([date ' ' mouse '- Reward (black), Omit (red)'])
            %savefig(fullfile(lg_out,img_fn, [img_fn '_repsByCrus.fig']))
        end
    end
end
    
        [x y] = ginput;
        
        
        xvals = 1:sz(2);
        yvals = 1:sz(1);
        
        y(1) = 1;
        y(end) = sz(1);
        x_int = interp1(y,x,yvals);
        
        split_img = zeros(size(img_ref));
        for i = yvals
            split_img(i,round(x_int(i)):sz(2)) = 1;
        end
        figure; imagesc(split_img)
        
        if exist(fullfile(lg_out,img_fn, [img_fn '_ROI_TCs.mat']))
            load(fullfile(lg_out,img_fn, [img_fn '_ROI_TCs.mat']))
        elseif exist(fullfile(jake_dir,img_fn, [img_fn '_ROI_TCs.mat']))
            load(fullfile(jake_dir,img_fn, [img_fn '_ROI_TCs.mat']))
        end
        nmask = size(mask3D,3);
        [maskCat maskCat_map] = splitMasks(split_img, mask3D, nmask);
        figure; imagesc(maskCat_map)

        save(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']), 'split_img', 'x', 'y','maskCat','maskCat_map');
    end
end