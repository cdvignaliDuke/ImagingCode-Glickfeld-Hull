
%all dF/F for bestSFTF
figure;
subplot(2,2,1)
col = strvcat('c', 'k', 'r');
x = [];
for iArea = 1:3
    [H, stats, xCDF, yCDF] = cdfplot_LG(squeeze(max(all_fits_dir(iArea).data_bestSFTF,[],1)));
    plot(xCDF, yCDF, col(iArea,:));
    hold on
    x= [x; size(all_fits_dir(iArea).x_bestSFTF(:,3),1)];
end
xlim([0 1.5])
legend(num2str(x));
xlabel('dF/F')
fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_dFoverF_bestSFTF.ps']);
print(gcf, '-depsc', fn_out);

%all dF/F at each SF/TF for each area
figure;
for iArea = 1:3
    subplot(2,3,iArea)
    [H, stats, xCDF, yCDF] = cdfplot_LG(max(all_fits_dir(iArea).data_bestSFTF_2Hz_pt16cpd,[],1));
    plot(xCDF, yCDF, 'b');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(max(all_fits_dir(iArea).data_bestSFTF_2Hz_pt04cpd,[],1));
    plot(xCDF, yCDF, 'c');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(max(all_fits_dir(iArea).data_bestSFTF_8Hz_pt16cpd,[],1));
    plot(xCDF, yCDF, 'g');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(max(all_fits_dir(iArea).data_bestSFTF_8Hz_pt04cpd,[],1));
    plot(xCDF, yCDF, 'r');
    hold on
    title(areas(iArea,:))
    xlabel('dF/F')
    xlim([0 1.5])
end
for iArea = 1:3
    subplot(2,3,3+iArea)
    [H, stats, xCDF, yCDF] = cdfplot_LG(max(all_fits_dir(iArea).data_2Hz_pt16cpd,[],1));
    plot(xCDF, yCDF, 'b');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(max(all_fits_dir(iArea).data_2Hz_pt04cpd,[],1));
    plot(xCDF, yCDF, 'c');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(max(all_fits_dir(iArea).data_8Hz_pt16cpd,[],1));
    plot(xCDF, yCDF, 'g');
    hold on
    [H, stats, xCDF, yCDF] = cdfplot_LG(max(all_fits_dir(iArea).data_8Hz_pt04cpd,[],1));
    plot(xCDF, yCDF, 'r');
    hold on
    title(areas(iArea,:))
    xlabel('dF/F')
    xlim([0 1.5])
end


fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_dFoverF_iArea_bestSFTF_iSFTF.ps']);
print(gcf, '-depsc', fn_out);

figure;
for iArea = 1:3
    subplot(2,2,1)
    scatter(max(all_fits_dir(iArea).data_bestSFTF_2Hz_pt04cpd,[],1), all_fits_dir(iArea).x_bestSFTF_2Hz_pt04cpd(:,3),'.',col(iArea,:));
    hold on
    xlabel('dF/F')
    ylabel('OSI')
    xlim([0 2])
    title('2Hz 0.04cpd')
    subplot(2,2,2)
    scatter(max(all_fits_dir(iArea).data_bestSFTF_2Hz_pt16cpd,[],1), all_fits_dir(iArea).x_bestSFTF_2Hz_pt16cpd(:,3),'.',col(iArea,:));
    hold on
    xlabel('dF/F')
    ylabel('OSI')
    xlim([0 2])
    title('2Hz 0.16cpd')
    subplot(2,2,3)
    scatter(max(all_fits_dir(iArea).data_bestSFTF_8Hz_pt04cpd,[],1), all_fits_dir(iArea).x_bestSFTF_8Hz_pt04cpd(:,3),'.',col(iArea,:));
    hold on
    xlabel('dF/F')
    ylabel('OSI')
    xlim([0 2])
    title('8Hz 0.04cpd')
    subplot(2,2,4)
    scatter(max(all_fits_dir(iArea).data_bestSFTF_8Hz_pt16cpd,[],1), all_fits_dir(iArea).x_bestSFTF_8Hz_pt16cpd(:,3),'.',col(iArea,:));
    hold on
    xlabel('dF/F')
    ylabel('OSI')
    xlim([0 2])
    title('8Hz 0.16cpd')
end

fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_dFoverF_OSI_bestSFTF_scatter_iSFTF.ps']);
print(gcf, '-depsc', fn_out);

figure;
for iArea = 1:3
    subplot(2,2,1)
    scatter(max(all_fits_dir(iArea).data_orituned_bestSFTF_2Hz_pt04cpd,[],1), all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt04cpd(:,8),'.',col(iArea,:));
    hold on
    xlabel('dF/F')
    ylabel('DSI')
    xlim([0 2])
    title('2Hz 0.04cpd')
    subplot(2,2,2)
    scatter(max(all_fits_dir(iArea).data_orituned_bestSFTF_2Hz_pt16cpd,[],1), all_fits_dir(iArea).x_orituned_bestSFTF_2Hz_pt16cpd(:,8),'.',col(iArea,:));
    hold on
    xlabel('dF/F')
    ylabel('DSI')
    xlim([0 2])
    title('2Hz 0.16cpd')
    subplot(2,2,3)
    scatter(max(all_fits_dir(iArea).data_orituned_bestSFTF_8Hz_pt04cpd,[],1), all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt04cpd(:,8),'.',col(iArea,:));
    hold on
    xlabel('dF/F')
    ylabel('DSI')
    xlim([0 2])
    title('8Hz 0.04cpd')
    subplot(2,2,4)
    scatter(max(all_fits_dir(iArea).data_orituned_bestSFTF_8Hz_pt16cpd,[],1), all_fits_dir(iArea).x_orituned_bestSFTF_8Hz_pt16cpd(:,8),'.',col(iArea,:));
    hold on
    xlabel('dF/F')
    ylabel('DSI')
    xlim([0 2])
    title('8Hz 0.16cpd')
end

fn_out = fullfile('\\zmey\storlab\users\Lindsey\Projects\HVAs\2P Axon Imaging', ['DIR_' num2str(P) 'P_' inj], ['DIR_' num2str(P) 'P_' inj '_dFoverF_DSI_bestSFTF_scatter_iSFTF.ps']);
print(gcf, '-depsc', fn_out);

figure;
for iArea = 1:3
    nexp = all_fits_dir(iArea).nexp;
    for iexp = 1:nexp
        n_pix = all_fits(iArea).expt(iexp).n(1,1);
        if size(all_fits_dir(iArea).expt(iexp).sfsf_tftf_grid,2)>1
            if find(all_fits_dir(iArea).expt(iexp).sfsf_tftf_grid(2,:) == 2)
                x = find(all_fits_dir(iArea).expt(iexp).sfsf_tftf_grid(2,:) == 2);
                if find(all_fits_dir(iArea).expt(iexp).sfsf_tftf_grid(1,:) == 0.04)
                    y = find(all_fits_dir(iArea).expt(iexp).sfsf_tftf_grid(1,:) == 0.04);
                    ind_A = intersect(x,y);
                    if find(all_fits_dir(iArea).expt(iexp).sfsf_tftf_grid(1,:) == 0.16)
                        z = find(all_fits_dir(iArea).expt(iexp).sfsf_tftf_grid(1,:) == 0.16);
                        ind_B = intersect(x,z);
                        for iCell = 1:n_pix
                            if all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit(ind_A) == 1
                                if all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit(ind_B) == 1
                                    subplot(2,3,iArea)
                                    if max(all_fits_dir(iArea).expt(iexp).bouton(iCell).dFoverF(:,ind_A),[],1)./max(all_fits_dir(iArea).expt(iexp).bouton(iCell).dFoverF(:,ind_B),[],1)>1
                                    scatter(all_fits_dir(iArea).expt(iexp).bouton(iCell).OSI(ind_A,:),all_fits_dir(iArea).expt(iexp).bouton(iCell).OSI(ind_B,:),'.','k');
                                    hold on
                                    else
                                    scatter(all_fits_dir(iArea).expt(iexp).bouton(iCell).OSI(ind_A,:),all_fits_dir(iArea).expt(iexp).bouton(iCell).OSI(ind_B,:),'.','r');
                                    hold on
                                    end
                                end
                            end
                            if all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit_orituned(ind_A) == 1
                                if all_fits_dir(iArea).expt(iexp).bouton(iCell).goodfit_orituned(ind_B) == 1
                                    subplot(2,3,iArea+3)
                                    if max(all_fits_dir(iArea).expt(iexp).bouton(iCell).dFoverF(:,ind_A),[],1)./max(all_fits_dir(iArea).expt(iexp).bouton(iCell).dFoverF(:,ind_B),[],1)>1
                                    scatter(all_fits_dir(iArea).expt(iexp).bouton(iCell).DSI(ind_A,:),all_fits_dir(iArea).expt(iexp).bouton(iCell).DSI(ind_B,:),'.','k');
                                    hold on
                                    else
                                    scatter(all_fits_dir(iArea).expt(iexp).bouton(iCell).DSI(ind_A,:),all_fits_dir(iArea).expt(iexp).bouton(iCell).DSI(ind_B,:),'.','r');
                                    hold on
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    subplot(2,3,iArea)
    title(areas(iArea,:))
    xlabel('OSI 0.04 cpd')
    ylabel('OSI 0.16 cpd')
    axis square
    subplot(2,3,iArea+3)
    title(areas(iArea,:))
    xlabel('DSI 0.04 cpd')
    ylabel('DSI 0.16 cpd')
    axis square
end

