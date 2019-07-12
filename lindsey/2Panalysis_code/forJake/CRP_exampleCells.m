CRP_expt_list_all
area_list = {'C1','C2','LS'};
for id = 1:3
    lg_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lindsey\Analysis\2P\Jake';
    share_out = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\public\ClassicConditioningPaper';
    nexp = size(expt(id).date,1);
    for iexp = 1:nexp

        mouse = expt(id).mouse(iexp,:);
        date = expt(id).date(iexp,:);
        run = expt(id).run(iexp,:);
        fprintf([date ' ' mouse '\n'])
        img_fn = [date '_' mouse];

        load(fullfile(lg_out,img_fn, [img_fn '_targetAlign.mat']))

        if exist(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            load(fullfile(lg_out,img_fn, [img_fn '_splitImage.mat']))
            area_temp = zeros(1,size(targetAligndFoverF,2));
            for i = 1:2
                if ~isnan(cell2mat(expt(id).areas{iexp}(i)))
                    x = find(strcmp(area_list, expt(id).areas{iexp}(i)));
                    ind = find(maskCat==i);
                    area_temp(1,ind) = x;
                else
                    ind = find(maskCat==i);
                    area_temp(1,ind) = NaN;
                end
            end
        else
            x = find(strcmp(area_list, expt(id).areas{iexp}));
            area_temp = ones(1,size(targetAligndFoverF,2)).*x;
        end
        areas = unique(area_temp(~isnan(area_temp)));
        for ia = areas
            if ~exist([share_out '\CC_examples\Day' num2str(id) '\' area_list{ia}])
                mkdir([share_out '\CC_examples\Day' num2str(id) '\' area_list{ia}])
            end
            ind = find(area_temp==ia);
            n = length(ind)-32;
            if n<0
                continue
            elseif n>16
                n = 16;
            end
            ntrials = length(ind_rew);
            if ntrials>50
                ntrials =50;
            end
            figure;
            for iC = 1:n
                i = ind(iC+32);
                subplot(4,4,iC)
                plot(tt,squeeze(targetAligndFoverF(:,i,ind_rew(1:ntrials))),'c')
                hold on
                plot(tt,nanmean(targetAligndFoverF(:,i,ind_rew(1:ntrials)),3),'b')
                ylabel('DF/F')
                xlabel('Time (ms)')
                title(['Cell ' num2str(iC+32)])
            end
            suptitle([mouse ' ' date])
            savefig([share_out '\CC_examples\Day' num2str(id) '\' area_list{ia} '\' mouse '_' date '_exampleCells_33-48.fig'])
        end
    end
end
% 
% for iC = 1:10
%     figure;
%     for i = 1:16
%         subplot(4,4,i)
%         plot(tt,squeeze(targetAligndFoverF(:,iC,i+10)))
%         hold on
%         ylim([-0.5 1])
%         if sum(targetAlign_events(:,iC,i+10),1)>0
%             scatter(tt(find(targetAlign_events(:,iC,i+10))),.75.*ones(sum(targetAlign_events(:,iC,i+10),1),1),'ok')
%         end
%         hold on
%         if sum(lickCueAlign(:,ind_rew(i+10)),1)>0
%             scatter(tt(find(lickCueAlign(:,i+10))),-0.25.*ones(sum(lickCueAlign(:,i+10),1),1),'or')
%         end
%         if find(ind_rew==i+10)
%             title('Rew')
%         else
%             title('Omit')
%         end
%     end
%     suptitle([mouse ' ' date '- Cell: #' num2str(iC) '- trial 11:26'])
%     savefig([share_out '\CC_examples\Day' num2str(id) '_' mouse '_' date '_exampleCell' num2str(iC) '.fig'])
% end
% 
% 
% 
