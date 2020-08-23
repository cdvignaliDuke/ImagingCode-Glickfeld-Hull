clc; clear all; close all;
doRedChannel = 0;
ds = 'CrossOriRandPhase_ExptList';
eval(ds)
rc = behavConstsAV;
frame_rate = 30;
nexp = size(expt,2);
phase_corr_all = zeros(8,2,nexp);

%%
for iexp = 1:nexp
    mouse = expt(iexp).mouse;
    date = expt(iexp).date;
    area = expt(iexp).img_loc{1};
    ImgFolder = expt(iexp).dirFolder;
    time = expt(iexp).dirTime;
    nrun = length(ImgFolder);
    dir_run_str = catRunName(cell2mat(ImgFolder), nrun);

    LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

    fprintf([mouse ' ' date '\n'])

    %% load data

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_TCs.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_input.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_stimData.mat']))
    
    nOn = input.nScansOn;
    nOff = input.nScansOff;
    nCells = size(npSub_tc,2);
    nTrials = length(dir_mat)-1;
    dir_mat = dir_mat(:,1:end-1);
    tc_crop = npSub_tc(nOff-frame_rate+1:end-nOn-frame_rate,:);
    tc_trial = reshape(tc_crop,[nOn+nOff nTrials nCells]);
    tc_f = mean(tc_trial(1:frame_rate,:,:),1);
    tc_df = tc_trial-tc_f;
    tc_dfof = tc_df./tc_f;
    
    tc_dir_dfof = zeros(nOn+nOff,nCells,nDir);
    for idir = 1:nDir
        ind = find(dir_mat == dirs(idir));
        tc_dir_dfof(:,:,idir) = squeeze(mean(tc_dfof(:,ind,:),2));
    end
    
    tc_dir_dfof_long = reshape(permute(tc_dir_dfof,[1 3 2]),[nDir.*(nOn+nOff) nCells]);
    [pca_coeff pca_score pca_latent pca_tsq pca_exp] = pca(tc_dir_dfof_long);
    pca_dir_dfof = reshape(pca_score, [nOn+nOff nDir nCells]);
    gray_lut = hsv(nDir);
    figure;
    for idir = 1:nDir
        plot3(pca_dir_dfof(:,idir,1), pca_dir_dfof(:,idir,2), pca_dir_dfof(:,idir,3),'-ok', 'MarkerFaceColor',gray_lut(idir,:))
        hold on
    end
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_dirResp_PCAspace_3D.pdf']),'-dpdf','-fillpage')

    pca_dir = reshape(permute(pca_dir_dfof,[1 3 2]),[(nOn+nOff).*nCells nDir]);
    corr_dir = corrcoef(pca_dir);
    figure; imagesc(corr_dir);
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' dir_run_str], [date '_' mouse '_' dir_run_str '_dirResp_PCAspace_corr.pdf']),'-dpdf','-fillpage')

    ImgFolder = expt(iexp).coFolder;
    time = expt(iexp).coTime;
    nrun = length(ImgFolder);
    co_run_str = catRunName(cell2mat(ImgFolder), nrun);

    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_respData.mat']))
    load(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_stimData.mat']))
    
    test_alone = data_dfof_con_ph_tc_avg(:,:,1,2,1);
    mask_alone = data_dfof_con_ph_tc_avg(:,:,2,1,1);
    all_plaid = squeeze(data_dfof_con_ph_tc_avg(:,:,2,2,:));
    
    nF = size(test_alone,1);
    test_pca = test_alone*pca_coeff;
    mask_pca = mask_alone*pca_coeff;
    plaid_pca = zeros(nF,nCells,nMaskPhas);
    for ip = 1:nMaskPhas
        plaid_pca(:,:,ip) = all_plaid(:,:,ip)*pca_coeff;
    end
    
%     figure;
%     for ip = 1:nMaskPhas
%         subplot(3,3,ip)
%         plot3(test_pca(:,1),test_pca(:,2), test_pca(:,3),'-o')
%         hold on
%         plot3(mask_pca(:,1),mask_pca(:,2),mask_pca(:,3),'-o')
%         plot3(plaid_pca(:,1,ip),plaid_pca(:,2,ip),plaid_pca(:,3,ip),'-o')
%     end
%     subplot(3,3,9)
%     plot3(test_pca(:,1),test_pca(:,2), test_pca(:,3),'-o')
%     hold on
%     plot3(mask_pca(:,1),mask_pca(:,2),mask_pca(:,3),'-o')
%     for ip = 1:2:nMaskPhas
%         plot3(plaid_pca(:,1,ip),plaid_pca(:,2,ip),plaid_pca(:,3,ip),'-o')
%     end
%     suptitle([mouse ' ' date])
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_phaseResp_PCAspace_3D.pdf']),'-dpdf','-fillpage')
% 
%     figure;
%     for ip = 1:nMaskPhas
%         subplot(3,3,ip)
%         plot(test_pca(:,1),test_pca(:,2),'-o')
%         hold on
%         plot(mask_pca(:,1),mask_pca(:,2),'-o')
%         plot(plaid_pca(:,1,ip),plaid_pca(:,2,ip),'-o')
%     end
%     subplot(3,3,9)
%     plot(test_pca(:,1),test_pca(:,2),'-o')
%     hold on
%     plot(mask_pca(:,1),mask_pca(:,2),'-o')
%     for ip = 1:2:nMaskPhas
%         plot(plaid_pca(:,1,ip),plaid_pca(:,2,ip),'-o')
%     end
%     suptitle([mouse ' ' date])
%     print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_phaseResp_PCAspace_2D.pdf']),'-dpdf','-fillpage')
    
    phase_corr = zeros(nMaskPhas,2);
    sz = size(test_pca);
    test_pca_long = reshape(test_pca,[sz(1).*nCells 1]);
    mask_pca_long = reshape(mask_pca,[sz(1).*nCells 1]);
    plaid_pca_long = reshape(plaid_pca,[sz(1).*nCells nMaskPhas]);
    for ip = 1:nMaskPhas
        phase_corr(ip,1) = triu2vec(corrcoef(test_pca_long, plaid_pca_long(:,ip)));
        phase_corr(ip,2) = triu2vec(corrcoef(mask_pca_long, plaid_pca_long(:,ip)));
    end
    phase_corr_all(:,:,iexp) = phase_corr;
    figure;
    plot(maskPhas,phase_corr(:,1))
    hold on
    plot(maskPhas,phase_corr(:,2))
    ylabel('Correlation')
    xlabel('Phase (deg)')
    legend({'test','mask'},'location','southeast')
    ylim([0 1])
    title([mouse ' ' date])
    print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_test&maskCorr.pdf']),'-dpdf','-fillpage')
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' co_run_str], [date '_' mouse '_' co_run_str '_pcaResp.mat']),'test_pca','mask_pca','plaid_pca','pca_score','pca_coeff','pca_latent','pca_exp','phase_corr')
end
summaryDir = fullfile(LG_base, 'Analysis', '2P', 'CrossOri', 'RandPhaseSummary');

figure;
errorbar(maskPhas,mean(phase_corr_all(:,1,:),3),std(phase_corr_all(:,1,:),[],3)./sqrt(nexp),'-o')
hold on
errorbar(maskPhas,mean(phase_corr_all(:,2,:),3),std(phase_corr_all(:,2,:),[],3)./sqrt(nexp),'-o')
ylabel('Correlation')
xlabel('Phase (deg)')
legend({'test','mask'},'location','southeast')
ylim([0 1])

print(fullfile(summaryDir, 'test&maskCorr_PCA.pdf'),'-dpdf','-fillpage')
    