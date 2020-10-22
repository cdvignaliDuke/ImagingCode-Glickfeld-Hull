clearvars
%close all
CRP_expt_list_mike;
pairings = Load_CRP_List();

for  j = 1:size(pairings,2)
id = pairings{1,j};
exp_subset = pairings{2,j};
nexp = size(expt(id).date,1);

mike_out = 'H:\home\mike\Analysis\CC 2P';
behavior_dir = 'H:\home\mike\Analysis\CRP Figures\Sessions & Summaries';


fprintf(['Day ' num2str(id) '\n'])
    for iexp = exp_subset
        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];

        if exist(fullfile(mike_out,img_fn, [img_fn '_reg.mat']))
            load(fullfile(mike_out,img_fn, [img_fn '_reg.mat']))
        end
            
        if exist(fullfile(mike_out,img_fn, [img_fn '_ROI_TCs.mat']))
            load(fullfile(mike_out,img_fn, [img_fn '_ROI_TCs.mat']))
            if exist(fullfile(mike_out, img_fn, ['AVG_img' mouse '.jpg']))
                figure; imshow(fullfile(mike_out, img_fn, ['AVG_img' mouse '.jpg'])); hold on;
                for i  = 1:size(mask3D,3)
                    bound = cell2mat(bwboundaries(mask3D(:,:,i)));
                    randcolor = rand(1,4);
                    plot(bound(:,2),bound(:,1),'.','color',randcolor); hold on;
                    text(mean(bound(:,2)),mean(bound(:,1)), ...
                        num2str(i), 'color', 'y', 'FontSize', 8);
                end
            end
        end
        
        if exist(fullfile(mike_out,img_fn, [img_fn '_splitImage.mat']))
            load(fullfile(mike_out,img_fn, [img_fn '_splitImage.mat']))
        else 
            sz = size(img_ref);
            figure;
            imagesc(img_ref)

            [x,y] = ginput;


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
            nmask = size(mask3D,3);
            [maskCat, maskCat_map] = splitMasks(split_img, mask3D, nmask);
            figure; imagesc(maskCat_map)
            save(fullfile(mike_out,img_fn, [img_fn '_splitImage.mat']), 'split_img', 'x', 'y','maskCat','maskCat_map');
        end
        
        load(fullfile(mike_out,img_fn, [img_fn '_targetAlign.mat']))
        load(fullfile(behavior_dir,['img' mouse '_sessions'], [date '_img' mouse '_various_licking_behaviors_by_trial.mat']))
        block2_trials_with_predictive_licking = intersect(trials_where_licking_preceded_reward,ind_block2); 
        ind_btwpl = block2_trials_with_predictive_licking;
        indL = find(maskCat==1);
        indR = find(maskCat==2);
        
        figure;
        shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,:,~ind_btwpl),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,:,~ind_btwpl),3),[],2)./sqrt(length(indL)+length(indR))).*(1000./frameRateHz),'r');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['CS- Licking Before Reward Time n=' num2str(length(indL)+length(indR))])
            vline(770)
            
            
 %%           
            figure;
            subplot(2,2,1)
            %shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,indL,ind_rew),3),2), (nanstd(nanmean(targetAligndFoverF(:,indL,ind_rew),3),[],2))./sqrt(length(indL)),'k');
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indL,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indL,ind_rew),3),[],2)./sqrt(length(indL))).*(1000./frameRateHz),'k');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['CS+ Left side- n=' num2str(length(indL))])
            vline(650)
            subplot(2,2,2)
            %shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,indR,ind_rew),3),2), (nanstd(nanmean(targetAligndFoverF(:,indR,ind_rew),3),[],2))./sqrt(length(indR)),'k');
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indR,ind_rew),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indR,ind_rew),3),[],2)./sqrt(length(indR))).*(1000./frameRateHz),'k');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            vline(770)
            title(['CS+ Right side- n=' num2str(length(indR))])
            subplot(2,2,3)
            %shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,indL,ind_block2),3),2), (nanstd(nanmean(targetAligndFoverF(:,indL,ind_block2),3),[],2))./sqrt(length(indL)),'r');
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indL,ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indL,ind_block2),3),[],2)./sqrt(length(indL))).*(1000./frameRateHz),'r');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['CS- Left side- n=' num2str(length(indL))])
            vline(770)
            subplot(2,2,4)
            %shadedErrorBar(tt,nanmean(nanmean(targetAligndFoverF(:,indR,ind_block2),3),2), (nanstd(nanmean(targetAligndFoverF(:,indR,ind_block2),3),[],2))./sqrt(length(indR)),'r');
            shadedErrorBar(tt,nanmean(nanmean(targetAlign_events(:,indR,ind_block2),3),2).*(1000./frameRateHz), (nanstd(nanmean(targetAlign_events(:,indR,ind_block2),3),[],2)./sqrt(length(indR))).*(1000./frameRateHz),'r');
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['CS- Right side- n=' num2str(length(indR))])
            vline(770)
            suptitle([date ' ' mouse '- Reward (black), None (red)'])
            
                       
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            figure;

            %subplot(2,2,1)
            for i = 1:length(indL)
                cell_mean = nanmean(targetAlign_events(:,indL(i),ind_rew),3).*(1000./frameRateHz); 
                p = plot(tt,cell_mean); hold on
                color = get(p, 'Color');
                text(tt(1)-200,cell_mean(1),num2str(indL(i)),'FontSize',6.,'color',color);
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['CS+ Left side- n=' num2str(length(indL))])
            vline(770)
            
            figure; %subplot(2,2,2)
            for i = 1:length(indR)
                cell_mean = nanmean(targetAlign_events(:,indR(i),ind_rew),3).*(1000./frameRateHz); 
                p = plot(tt,cell_mean); hold on
                color = get(p, 'Color');
                text(tt(1)-200,cell_mean(1),num2str(indR(i)),'FontSize',6.,'color',color);
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            vline(650)
            title(['CS+ Right side- n=' num2str(length(indR))])
            
            figure; %subplot(2,2,3)
            for i = 1:length(indL)
                cell_mean = nanmean(targetAlign_events(:,indL(i),ind_block2),3).*(1000./frameRateHz); 
                p = plot(tt,cell_mean); hold on
                color = get(p, 'Color');
                text(tt(1)-200,cell_mean(1),num2str(indL(i)),'FontSize',6.,'color',color);
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['CS- Left side- n=' num2str(length(indL))])
            vline(770)
            
            figure; %subplot(2,2,4)
            for i = 1:length(indR)
                cell_mean = nanmean(targetAlign_events(:,indR(i),ind_block2),3).*(1000./frameRateHz); 
                p = plot(tt,cell_mean); hold on
                color = get(p, 'Color');
                text(tt(1)-200,cell_mean(1),num2str(indR(i)),'FontSize',6.,'color',color);
                
            end
            xlabel('Time from cue')
            ylabel('Spike rate (Hz)')
            title(['CS- Right side- n=' num2str(length(indR))])
            vline(770)
            %suptitle([date ' ' mouse 'Cells Averaged Over Trials'])
            %savefig(fullfile(lg_out,img_fn, [img_fn '_repsByCrus.fig']))
        %end
    end
end

