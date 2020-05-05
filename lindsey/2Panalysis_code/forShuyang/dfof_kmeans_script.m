%%dfof_kmeans_script
clear;
days = {'191220_img1042'};
basedir = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff';
datadir = fullfile(basedir, 'home\shuyang\Analysis\2P_MovingDots_Analysis\imaging_analysis');
%set output directory
outbase = fullfile(basedir, 'public\MetaKMeans2020');
%timecourse file name base
tcstr = '_000_nPCA200_mu0.3_nIC200_thresh97_TCave.mat';
%mask file name base
maskstr = '_000_nPCA200_mu0.3_nIC200_thresh97_coor0.8_mask3D_final.mat';
for iexp = 1
    load(fullfile(datadir, days{iexp}, 'getTC', [days{iexp} tcstr]));
    load(fullfile(datadir, days{iexp}, 'getTC', [days{iexp} maskstr]));
    if ~exist(fullfile(outbase, days{iexp}))
        mkdir(fullfile(outbase, days{iexp}))
    end
    
    %transform TC to dfof
    meanF = mean(tc_avg, 1);
    df_f = bsxfun(@minus, tc_avg, meanF);
    df_f = bsxfun(@rdivide, df_f, meanF);
    max_dff = max(df_f, [], 1);
    df_f_norm = bsxfun(@rdivide, df_f, max_dff);
    
    
    mask_cell = zeros(size(mask3D(:,:,1)));
    for i = 1:size(tc_avg,2)
        mask_cell(find(mask3D(:,:,i))) = i;
    end
    
    % no seeds
    sc = 0;
%     % code for making seeds
%     h = figure;
%     imagesc(mask_cell)
%     [i j] = getpts(h);
%     nsc = size(i,1);
%     scs = zeros(1,nsc);
%     sc_mask = zeros(size(mask_cell));
%     for isc = 1:nsc
%         scs(isc)= mask_cell(round(j(isc)),round(i(isc)));
%         sc_mask(find(mask_cell==scs(isc))) = isc;
%     end
%     sc = df_f_norm(:,scs)';
    

    %dendrite correlation matrix
    nDend = size(tc_avg,2);
    r_dend = zeros(nDend);
    for id = 1:nDend
        for id2 = 1:nDend
            r_dend(id,id2) = triu2vec(corrcoef(tc_avg(:,id),tc_avg(:,id2)));
        end
    end
    
    figure;
    for i = 2:7
    idx = kmeans(r_dend,17);
    [x idx_sort] = sort(idx,'ascend');
    ind = find(diff(x));
    r_dend_sort = r_dend(idx_sort,idx_sort);
    subplot(3,2,i-1)
    imagesc(r_dend_sort)
    hline(ind+.5)
    vline(ind+.5)
    end
    
    print(fullfile(outbase, days{iexp}, [days{iexp} '_sortMatrix.pdf']),'-dpdf', '-bestfit')
    
    % run meta k means
    [allclusters, allclusters_init, centroidcorr, dendmem, dunnsinitial]=meta_k_means_LG(df_f_norm', 'correlation', sc);
    
    
    %plot output of meta k means
    figure
    mask_cluster = zeros(size(mask_cell));
    [n n2] = subplotn(size(allclusters,1)+2);
    for i = 1:size(allclusters,1)
        mask_cluster_temp = zeros(size(mask_cell));
        for ii = 1:size(allclusters{i,1},2)
            mask_cluster(find(mask_cell ==  allclusters{i,1}(ii))) = i;
            mask_cluster_temp(find(mask_cell ==  allclusters{i,1}(ii))) = 1;
        end
        subplot(n,n2,i)
        imagesc(mask_cluster_temp)
    end
    
    subplot(n,n2,i+1); imagesc(mask_cluster)
    
    % if using seeds
    if size(scs,2)>1
        subplot(n,n2,i+2); imagesc(sc_mask)
    end
    print(fullfile(outbase, days{iexp}, [days{iexp} '_clusters_sc.pdf']),'-dpdf')
    
    c_order = [];
    b = 0.5;
    for ic = 1:size(allclusters,1)
        c_order = [c_order allclusters{ic,1}];
        b = [b b(end)+length(allclusters{ic,1})];
    end
    figure; imagesc(r_dend(c_order,c_order))
    hline(b)
    vline(b)
    
    
    save(fullfile(outbase, days{iexp}, [days{iexp} '_clusters_sc.mat']), 'allclusters', 'centroidcorr', 'dendmem', 'dunnsinitial', 'mask_cell', 'mask_cluster', 'sc');
end
