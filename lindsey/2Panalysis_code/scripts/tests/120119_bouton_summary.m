fn_out = fullfile(sum_base, [matrix '_' num2str(P) 'P_' inj '_axon_summary'], 'all_fits.mat');
load(all_fits);
        
areas = ['PM'; 'LM'; 'AL'];
for iArea = 1:3;
    P = 2;
    matrix = 'SF5xTF5';
    inj = 'V1';
    image = areas(iArea,:);

    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    
    for iexp = 1:nexp;

        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};
        zoom = exp_list.zoom_mat{iexp};
        
        base = 'G:\users\lindsey\analysisLG\active mice';
        outDir = fullfile(base, mouse, date);

        n_pix = all_fits(iArea).expt(iexp).n(1);

        fn_stim = fullfile(outDir, 'analysis', [date '_' mouse '_run' num2str(userun) '_allstim.mat']);
        load(fn_stim);

        
        siz = size(stim_off);
        sq_pix_off = zeros(9,9,n_pix);
        sq_pix_on = zeros(9,9,n_pix);
        good_ind = [];
        bad_ind = [];
        highdF_ind = [];
        lowdF_ind = [];
        for iCell = 1:n_pix
            pos = all_fits(iArea).expt(iexp).bouton(iCell).pos;
            suby = [pos(1)-4:pos(1)+4];
            subx = [pos(2)-4:pos(2)+4];
            sq_pix_off(:,:,iCell) = mean(stim_off(suby,subx,:),3);
            sq_pix_on(:,:,iCell) = mean(stim_on(suby,subx,:),3);
            if all_fits(iArea).expt(iexp).bouton(iCell).goodfit == 1
                good_ind = [good_ind iCell];
                if all_fits(iArea).expt(iexp).bouton(iCell).dF_fit>.3
                    highdF_ind = [highdF_ind iCell];
                else
                    lowdF_ind = [lowdF_ind iCell];
                end
            elseif all_fits(iArea).expt(iexp).bouton(iCell).goodfit < 1
                bad_ind = [bad_ind iCell];
            end
        end
        bouton_off = mean(sq_pix_off,3);
        bouton_on = mean(sq_pix_on,3);
        bouton_off_good = mean(sq_pix_off(:,:,good_ind),3);
        bouton_off_bad = mean(sq_pix_off(:,:,bad_ind),3);
        bouton_off_highdF = mean(sq_pix_off(:,:,highdF_ind),3);
        bouton_off_lowdF = mean(sq_pix_off(:,:,lowdF_ind),3);
        
        fn_boutons = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_boutons.mat']);
        save(fn_boutons, 'sq_pix_off', 'sq_pix_on', 'bouton_off', 'bouton_on', 'bouton_off_good', 'bouton_off_bad', 'bouton_off_highdF', 'bouton_off_lowdF', 'good_ind', 'bad_ind', 'lowdF_ind', 'highdF_ind');
        
        figure;
        max_F = max(max(bouton_on,[],2),[],1);
        subplot(3,2,1); imagesq(bouton_off); caxis([0 max_F]); title(['off (' num2str(n_pix) ')']);
        subplot(3,2,2); imagesq(bouton_on); caxis([0 max_F]); title(['on (' num2str(n_pix) ')']);
        subplot(3,2,3); imagesq(bouton_off_good); caxis([0 max_F]); title(['good fit (' num2str(length(good_ind)) ')']);
        subplot(3,2,4); imagesq(bouton_off_bad); caxis([0 max_F]); title(['bad fit (' num2str(length(bad_ind)) ')']);
        subplot(3,2,5); imagesq(bouton_off_highdF); caxis([0 max_F]); title(['dF/F>0.3(' num2str(length(highdF_ind)) ')']);
        subplot(3,2,6); imagesq(bouton_off_lowdF); caxis([0 max_F]);title(['dF/F<0.3 (' num2str(length(lowdF_ind)) ')']);
        colormap(gray)
        
        fn_boutonfig = fullfile(anal_base, areas(iArea,:), [date '_' mouse '_run' num2str(userun) '_boutons.pdf']);
         print(gcf, '-dpdf', fn_out);
    end
end

for iArea = 1:3
    P = 2;
    matrix = 'SF5xTF5';
    inj = 'V1';
    image = areas(iArea,:);

    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);
    z1_pix_off = [];
    z1_pix_on = [];
    z1pt5_pix_off = [];
    z1pt5_pix_on = [];
    z1_pix_off_good = [];
    z1_pix_off_bad = [];
    z1_pix_off_highdF = [];
    z1_pix_off_lowdF = [];
    z1pt5_pix_off_good = [];
    z1pt5_pix_off_bad = [];
    z1pt5_pix_off_highdF = [];
    z1pt5_pix_off_lowdF = [];
    for iexp = 1:nexp;
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};
        zoom = exp_list.zoom_mat{iexp};
        
        base = 'G:\users\lindsey\analysisLG\active mice';
        outDir = fullfile(base, mouse, date);
        fn_boutons = fullfile(outDir,'analysis',[date '_' mouse '_run' num2str(userun) '_boutons.mat']);
        load(fn_boutons);
        if zoom == 1;
            z1_pix_off = cat(3,z1_pix_off,sq_pix_off);
            z1_pix_on = cat(3,z1_pix_on,sq_pix_on);
            z1_pix_off_good = cat(3,z1_pix_off_good,sq_pix_off(:,:,good_ind));
            z1_pix_off_bad = cat(3,z1_pix_off_bad,sq_pix_off(:,:,bad_ind));
            z1_pix_off_highdF = cat(3,z1_pix_off_highdF,sq_pix_off(:,:,highdF_ind));
            z1_pix_off_lowdF = cat(3,z1_pix_off_lowdF,sq_pix_off(:,:,lowdF_ind));
            
        elseif zoom ==1.5;
            z1pt5_pix_off = cat(3,z1pt5_pix_off,sq_pix_off);
            z1pt5_pix_on = cat(3,z1pt5_pix_on,sq_pix_on);
            z1pt5_pix_off_good = cat(3,z1pt5_pix_off_good,sq_pix_off(:,:,good_ind));
            z1pt5_pix_off_bad = cat(3,z1pt5_pix_off_bad,sq_pix_off(:,:,bad_ind));
            z1pt5_pix_off_highdF = cat(3,z1pt5_pix_off_highdF,sq_pix_off(:,:,highdF_ind));
            z1pt5_pix_off_lowdF = cat(3,z1pt5_pix_off_lowdF,sq_pix_off(:,:,lowdF_ind));
        end
    end

    z1_off = mean(z1_pix_off,3);
    z1_on = mean(z1_pix_on,3);
    z1pt5_off = mean(z1pt5_pix_off,3);
    z1pt5_on = mean(z1pt5_pix_on,3);
    z1_off_bad = mean(z1_pix_off_bad,3);
    z1_off_good = mean(z1_pix_off_good,3);
    z1_off_highdF = mean(z1_pix_off_highdF,3);
    z1_off_lowdF = mean(z1_pix_off_lowdF,3);
    z1pt5_off_bad = mean(z1pt5_pix_off_bad,3);
    z1pt5_off_good = mean(z1pt5_pix_off_good,3);
    z1pt5_off_highdF = mean(z1pt5_pix_off_highdF,3);
    z1pt5_off_lowdF = mean(z1pt5_pix_off_lowdF,3);
    
    
    figure;
    max_F = max(max(z1_on,[],2),[],1);
    subplot(3,2,1); imagesq(z1_off); caxis([0 max_F]); title(['off (' num2str(size(z1_pix_off,3)) ')']);
    subplot(3,2,2); imagesq(z1_on); caxis([0 max_F]); title(['on (' num2str(size(z1_pix_on,3)) ')']);
    subplot(3,2,3); imagesq(z1_off_good); caxis([0 max_F]); title(['good fit (' num2str(size(z1_pix_off_good,3)) ')']);
    subplot(3,2,4); imagesq(z1_off_bad); caxis([0 max_F]); title(['bad fit (' num2str(size(z1_pix_off_bad,3)) ')']);
    subplot(3,2,5); imagesq(z1_off_highdF); caxis([0 max_F]); title(['dF/F>0.3 (' num2str(size(z1_pix_off_highdF,3)) ')']);
    subplot(3,2,6); imagesq(z1_off_lowdF); caxis([0 max_F]);title(['dF/F<0.3 (' num2str(size(z1_pix_off_lowdF,3)) ')']);
    colormap(gray)
    suptitle(['All boutons ' areas(iArea,:) ' z=1' ]);
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' areas(iArea) '_z1_boutons.pdf']);
    print(gcf, '-dpdf', fn_out);
    figure;
    max_F = max(max(z1pt5_on,[],2),[],1);
    subplot(3,2,1); imagesq(z1pt5_off); caxis([0 max_F]); title(['off (' num2str(size(z1pt5_pix_off,3)) ')']);
    subplot(3,2,2); imagesq(z1pt5_on); caxis([0 max_F]); title(['on (' num2str(size(z1pt5_pix_on,3)) ')']);
    subplot(3,2,3); imagesq(z1pt5_off_good); caxis([0 max_F]); title(['good fit (' num2str(size(z1pt5_pix_off_good,3)) ')']);
    subplot(3,2,4); imagesq(z1pt5_off_bad); caxis([0 max_F]); title(['bad fit (' num2str(size(z1pt5_pix_off_bad,3)) ')']);
    subplot(3,2,5); imagesq(z1pt5_off_highdF); caxis([0 max_F]); title(['dF/F>0.3 (' num2str(size(z1pt5_pix_off_highdF,3)) ')']);
    subplot(3,2,6); imagesq(z1pt5_off_lowdF); caxis([0 max_F]);title(['dF/F<0.3 (' num2str(size(z1pt5_pix_off_lowdF,3)) ')']);
    colormap(gray)
    suptitle(['All boutons ' areas(iArea,:) ' z=1.5' ]);
    fn_out = fullfile(anal_base, [matrix '_' num2str(P) 'P_' inj '_' areas(iArea) '_z1pt5_boutons.pdf']);
    print(gcf, '-dpdf', fn_out);
end