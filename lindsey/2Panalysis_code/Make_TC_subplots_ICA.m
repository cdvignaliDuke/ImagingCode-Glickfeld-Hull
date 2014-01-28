    figure;
    areas = size(roi_avg,2);
    if run ==1;
        TC_norm_avg = squeeze(TC_norm_avg(:,:,:,1));
        TC_norm_sem = squeeze(TC_norm_sem(:,:,:,1));
        resp_avg = squeeze(resp_avg(:,:,1));
    elseif areas == 1;
        TC_norm_avg = squeeze(TC_norm_avg(:,:,:));
        TC_norm_sem = squeeze(TC_norm_sem(:,:,:));
        resp_avg = squeeze(resp_avg(:,:));
    end
        
    for iCond = 1:nCond;
        for iArea = 1:areas;
            h = subplot(sqrt(nCond)+1,sqrt(nCond),iCond);
            if areas==1
            errorbar(1:(nON+nOFF), TC_norm_avg(:,iCond), TC_norm_sem(:,iCond), col(iArea));
            mxy = max(max((TC_norm_avg),[],2),[],1);
            mny = min(min((TC_norm_avg),[],2),[],1);
            else    
            errorbar(1:(nON+nOFF), TC_norm_avg(:,iArea,iCond), TC_norm_sem(:,iArea,iCond), col(iArea));
            mxy = max(max(max((TC_norm_avg),[],3),[],2),[],1);
            mny = min(min(min((TC_norm_avg),[],3),[],2),[],1);
            end
            set(h,'XTick',[0:(nON+nOFF)/4:(nON+nOFF)]);
            set(h,'XTickLabel',(0:(nON+nOFF)/4:(nON+nOFF)));
            ylim([mny (mxy + 0.1*mxy)]);
            xlim([0 (nON+nOFF)+2]);
            hold on;
        end 
    end
    if blanks==1;
        iCond = nCond+1;
        for iArea = 1:areas;
            h=subplot(sqrt(nCond)+1,sqrt(nCond),iCond);
            if areas==1
            errorbar(1:(nON+nOFF), TC_norm_avg(:,iCond), TC_norm_sem(:,iCond), col(iArea));
            mxy = max(max((TC_norm_avg),[],2),[],1);
            mny = min(min((TC_norm_avg),[],2),[],1);
            else    
            errorbar(1:(nON+nOFF), TC_norm_avg(:,iArea,iCond), TC_norm_sem(:,iArea,iCond), col(iArea));
            my = max(max(max((TC_norm_avg),[],3),[],2),[],1);
            end
            set(h,'XTick',[0:(nON+nOFF)/4:(nON+nOFF)]);
            set(h,'XTickLabel',(0:(nON+nOFF)/4:(nON+nOFF)));
            ylim([mny (mxy + 0.1*mxy)]);
            xlim([0 (nON+nOFF)+2]);
            hold on;
        end
    end
    subplot(sqrt(nCond)+1,sqrt(nCond),iCond+1)
    imstretch(sm(:,:,iCell),[.5 .99],1.5);
    text(.8,.1,num2str(iCell),'fontsize',12,'color','w','fontweight','bold','unit','norm');
    %legend(area_list)