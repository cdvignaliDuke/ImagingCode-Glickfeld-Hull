
    cells = size(roi_avg,2);
    if run ==1;
        TC_avg = squeeze(TC_norm_avg(:,:,:,1));
        TC_sem = squeeze(TC_norm_sem(:,:,:,1));
        resp_avg = squeeze(resp_avg(:,:,1));
    end
    for iCell = 1:5;      
            figure;
            for iCond = 1:nCond;
            h = subplot(sqrt(nCond)+1,sqrt(nCond),iCond); 
            errorbar(1:(nON+nOFF), TC_norm_avg(:,iCell,iCond), TC_norm_sem(:,iCell,iCond));
            my = max(max(max((TC_norm_avg),[],3),[],2),[],1);
            set(h,'XTick',[0:(nON+nOFF)/4:(nON+nOFF)]);
            set(h,'XTickLabel',(0:(nON+nOFF)/4:(nON+nOFF)));
            %ylim([-.5 (my + 0.1*my)]);
            ylim([0 my]);
            xlim([0 (nON+nOFF)+2]);
            hold on;
        end
        if blanks==1;
            iCond = nCond+1;
            h=subplot(sqrt(nCond)+1,sqrt(nCond),iCond);
            errorbar(1:(nON+nOFF), TC_norm_avg(:,iCell,iCond), TC_norm_sem(:,iCell,iCond));
            my = max(max(max((TC_norm_avg),[],3),[],2),[],1);
            set(h,'XTick',[0:(nON+nOFF)/4:(nON+nOFF)]);
            set(h,'XTickLabel',(0:(nON+nOFF)/4:(nON+nOFF)));
            %ylim([-.5 (my + 0.1*my)]);
            ylim([0 my]);
            xlim([0 (nON+nOFF)+2]);
            hold on;
        end
    end
        

