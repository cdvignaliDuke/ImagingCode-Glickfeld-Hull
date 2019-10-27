%code for matching individual neurons between day 1 and day N. Then
%comparing neural responses between days. 

clear all
close all
CRP_expt_list
nd = size(expt,2);
behav_dir = '\\crash.dhe.duke.edu\data\home\andrew\Behavior\Data';
jake_dir = '\\crash.dhe.duke.edu\data\home\jake\Analysis\Cue_reward_pairing_analysis\2P';
lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
nexp = size(expt(1).date,1);
nexp2 = size(expt(2).date,1);
for iexp = 1:nexp
    for id  = 1:2
        if id == 1
            figure;
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
            cd(fullfile('\\crash.dhe.duke.edu\data\home\jake\Data\2P_imaging', [date '_img' mouse_name],['img' mouse_name]))
            clear global
            load(['img' mouse_name '_000_' run '.mat'])
            load(fullfile(lg_out, img_fn, [img_fn '_reg.mat']));
            fName = [mouse_name '_000_000'];
            data = sbxread(['img' mouse_name '_000_' run],0,1000);
            data = squeeze(data);
            [out img_reg1] = stackRegister(data,img_ref);
            subplot(3,3,id)
            img_reg1_avg = mean(img_reg1,3);
            imagesc(img_reg1_avg)
            colormap gray
            cmin = min(img_reg1_avg(:));
            cmax = max(img_reg1_avg(:))/3;
            if (cmax./cmin)<1.5
                cmax = max(img_reg1_avg(:));
            end
            clim([cmin cmax])
            title(date)
            load(fullfile(lg_out, img_fn, [img_fn '_ROI_TCs.mat']));
            subplot(3,3,id+3)
            imagesc(mask_flat)
            title(date)
            mask1 = sum(mask_flat,3);
            mask1(find(mask1>0)) = 1;
            mask3D_1 = mask3D;
            nmask1 = size(mask3D_1,3);
            load(fullfile(lg_out, img_fn, [img_fn '_targetAlign.mat']));
            ind_rew1 = ind_rew;
            ind_omit1 = ind_omit;
            targetAlign_events1 = targetAlign_events;
            targetAligndFoverF1 = targetAligndFoverF;
        elseif id == 2
            for iexp2 = 1:nexp2
                if strcmp(expt(id).mouse(iexp2,:),mouse)
                    break
                end
            end
            mouse = expt(id).mouse(iexp2,:);
            date = expt(id).date(iexp2,:);
            run = expt(id).run(iexp2,:);
            fprintf([mouse ' ' date '\n'])
            img_fn = [date '_' mouse];
            cd(fullfile('\\crash.dhe.duke.edu\data\home\jake\Data\2P_imaging', [date '_img' mouse_name],['img' mouse_name]))
            clear global
            load(['img' mouse_name '_000_' run '.mat'])
            load(fullfile(lg_out, img_fn, [img_fn '_reg.mat']));
            data = sbxread(['img' mouse_name '_000_' run],0,1000);
            data = squeeze(data);
            [out img_reg2] = stackRegister(data,img_ref);
            subplot(3,3,id)
            img_reg2_avg = mean(img_reg2,3);
            imagesc(img_reg2_avg)
            colormap gray
            cmin = min(img_reg2_avg(:));
            cmax = max(img_reg2_avg(:))/3;
            if (cmax./cmin)<1.5
                cmax = max(img_reg2_avg(:));
            end
            clim([cmin cmax])
            title(date)
            load(fullfile(lg_out, img_fn, [img_fn '_ROI_TCs.mat']));
            subplot(3,3,id+3)
            imagesc(mask_flat)
            title(date)
            mask2 = sum(mask_flat,3);
            mask2(find(mask2>0)) = 1;
            [D2toD1_out img_reg3] = stackRegister(img_reg2_avg,img_reg1_avg);
            subplot(3,3,id+1)
            imagesc(img_reg3)
            colormap gray
            cmin = min(img_reg3(:));
            cmax = max(img_reg3(:))/3;
            if (cmax./cmin)<1.5
                cmax = max(img_reg3(:));
            end
            clim([cmin cmax])
            title('Reg')
            [out2 mask_reg] = stackRegister_MA(mask2,[],[],D2toD1_out);
            subplot(3,3,6)
            mask_rgb = cat(3,mask1,mask_reg);
            mask_b = zeros(size(mask1));
            mask_rgb = cat(3,mask_rgb,mask_b);
            imagesc(mask_rgb)
            title('Reg overlay')
            nmask2 = size(mask3D,3);
            [out mask3D_2] = stackRegister_MA(mask3D,[],[],repmat(D2toD1_out,[nmask2 1]));
            nmask2 = size(mask3D_2,3);
            ind_overlap = cell(1,nmask2);
            pct_overlap = cell(1,nmask2);
            n_overlap = zeros(1,nmask2);
            pct_pix_overlap = [];
            ind_overlap_alt = cell(1,nmask2);
            n_overlap_alt = zeros(1,nmask2);
            pct_pix_overlap_alt = [];
            overlap_img = zeros(size(mask_rgb));
            overlap_id = nan(1,nmask2);
            for i = 1:nmask2
                temp = mask3D_2(:,:,i);
                temp(find(temp>0.5)) = 1;
                temp(find(temp<0.5)) = 0;
                mask3D_2(:,:,i) = temp;
                for ii = 1:nmask1
                    if find(mask3D_2(:,:,i)+mask3D_1(:,:,ii) == 2)
                        ind_overlap{i} = [ind_overlap{i} ii];
                        n_overlap(1,i) = n_overlap(1,i)+1; 
                        npix = length(find(mask3D_2(:,:,i)));
                        pct_pix_overlap = [pct_pix_overlap length(find(mask3D_2(:,:,i)+mask3D_1(:,:,ii) == 2))./npix];
                        pct_overlap{i} = [pct_overlap{i} pct_pix_overlap];
                        if pct_pix_overlap(end)>0.6
                            npix1 = length(find(mask3D_1(:,:,ii)));
                            pct_pix_overlap2 = length(find(mask3D_2(:,:,i)+mask3D_1(:,:,ii) == 2))./npix1;
                            if pct_pix_overlap2>0.6
                                overlap_img(:,:,1) = overlap_img(:,:,1) + mask3D_1(:,:,ii);
                                overlap_img(:,:,2) = overlap_img(:,:,2) + mask3D_2(:,:,i);
                                overlap_id(1,i) = ii;
                            end
                        end
                    end
                end
            end
            subplot(3,3,7)
            hist(n_overlap)
            ylabel('Dendrites')
            xlabel('Number overlapping')
            subplot(3,3,8)
            hist(pct_pix_overlap)
            ylabel('Dendrites')
            xlabel('Percent of overlap')
            xlim([0 1])
            title('Correct')
            subplot(3,3,9)
            imagesc(overlap_img)
            load(fullfile(lg_out, img_fn, [img_fn '_targetAlign.mat']));
            ind_rew2 = ind_rew;
            ind_omit2 = ind_omit;
            targetAlign_events2 = targetAlign_events;
            targetAligndFoverF2 = targetAligndFoverF;
        end
    end
    suptitle(mouse)
    print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_crossday\' mouse '_refs.pdf'],'-bestfit','-dpdf')
    print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_crossday\' mouse '.eps'],'-depsc')
%     save(fullfile(lg_out, img_fn, [img_fn '_D2toD1_overlap.mat']), 'D2toD1_out', 'ind_overlap', 'pct_overlap','overlap_id')
%     if sum(~isnan(overlap_id))
%         nmatch = sum(~isnan(overlap_id));
%         ind = find(~isnan(overlap_id));
%         [n n2] = subplotn(nmatch);
%         figure;
%         for i = 1:nmatch
%             subplot(n,n2,i)
%             shadedErrorBar(tt,nanmean(targetAlign_events1(:,overlap_id(ind(i)),ind_rew1),3).*1000/frameRateHz, (nanstd(targetAlign_events1(:,overlap_id(ind(i)),ind_rew1),[],3).*1000/frameRateHz)./sqrt(length(ind_rew2)),'k');
%             hold on
%             shadedErrorBar(tt,nanmean(targetAlign_events2(:,ind(i),ind_rew2),3).*1000/frameRateHz, (nanstd(targetAlign_events2(:,ind(i),ind_rew2),[],3).*1000/frameRateHz)./sqrt(length(ind_rew2)),'b');
%             xlabel('Time from cue (ms)')
%             ylabel('Spike rate (Hz)')
%             ylim([-1 8])
%             xlim([-500 2000])
%             vline([600], 'b')
%         end
%         suptitle([mouse ' ' date ' rewarded trials- D1 (black) vs D2 (blue)'])
%         print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_crossday\' mouse '_matchedCellsD1vD2_spikes_reward.pdf'],'-bestfit','-dpdf')
%         
%         figure;
%         for i = 1:nmatch
%             subplot(n,n2,i)
%             shadedErrorBar(tt,nanmean(targetAligndFoverF1(:,overlap_id(ind(i)),ind_rew1),3), (nanstd(targetAligndFoverF1(:,overlap_id(ind(i)),ind_rew1),[],3))./sqrt(length(ind_rew2)),'k');
%             hold on
%             shadedErrorBar(tt,nanmean(targetAligndFoverF2(:,ind(i),ind_rew2),3), (nanstd(targetAligndFoverF2(:,ind(i),ind_rew2),[],3))./sqrt(length(ind_rew2)),'b');
%             xlabel('Time from cue (ms)')
%             ylabel('dF/F')
%             ylim([-0.1 1])
%             xlim([-500 2000])
%             vline([600], 'b')
%         end
%         suptitle([mouse ' ' date ' rewarded trials- D1 (black) vs D2 (blue)'])
%         print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_crossday\' mouse '_matchedCellsD1vD2_dfoverf_reward.pdf'],'-bestfit','-dpdf')
% 
%         figure;
%         for i = 1:nmatch
%             subplot(n,n2,i)
%             shadedErrorBar(tt,nanmean(targetAlign_events1(:,overlap_id(ind(i)),ind_omit1),3).*1000/frameRateHz, (nanstd(targetAlign_events1(:,overlap_id(ind(i)),ind_omit1),[],3).*1000/frameRateHz)./sqrt(length(ind_omit2)),'k');
%             hold on
%             shadedErrorBar(tt,nanmean(targetAlign_events2(:,ind(i),ind_omit2),3).*1000/frameRateHz, (nanstd(targetAlign_events2(:,ind(i),ind_omit2),[],3).*1000/frameRateHz)./sqrt(length(ind_omit2)),'b');
%             xlabel('Time from cue (ms)')
%             ylabel('Spike rate (Hz)')
%             ylim([-1 8])
%             xlim([-500 2000])
%             vline([600], 'b')
%         end
%         suptitle([mouse ' ' date ' omission trials- D1 (black) vs D2 (blue)'])
%         print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_crossday\' mouse '_matchedCellsD1vD2_spikes_omission.pdf'],'-bestfit','-dpdf')
%         
%         figure;
%         for i = 1:nmatch
%             subplot(n,n2,i)
%             shadedErrorBar(tt,nanmean(targetAligndFoverF1(:,overlap_id(ind(i)),ind_omit1),3), (nanstd(targetAligndFoverF1(:,overlap_id(ind(i)),ind_omit1),[],3))./sqrt(length(ind_omit2)),'k');
%             hold on
%             shadedErrorBar(tt,nanmean(targetAligndFoverF2(:,ind(i),ind_omit2),3), (nanstd(targetAligndFoverF2(:,ind(i),ind_omit2),[],3))./sqrt(length(ind_omit2)),'b');
%             xlabel('Time from cue (ms)')
%             ylabel('dF/F')
%             ylim([-0.1 1])
%             xlim([-500 2000])
%             vline([600], 'b')
%         end
%         suptitle([mouse ' ' date ' omission trials- D1 (black) vs D2 (blue)'])
%         print(['\\crash.dhe.duke.edu\data\public\ClassicConditioningPaper\CC_crossday\' mouse '_matchedCellsD1vD2_dfoverf_omission.pdf'],'-bestfit','-dpdf')
% 
%     end
    %close all
end
            