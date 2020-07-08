%%dfof_kmeans_script
clear;
sessions = {'200305_img1049'};
days = {'1049-200305_1'};
image_analysis_base    = 'Z:\Analysis\motorizedWheel_Analysis\airpuff\imaging_analysis\';

for iexp = 1
    image_analysis_dest = [image_analysis_base, sessions{iexp}, '\'];
    outbase = [image_analysis_base, sessions{iexp}, '\metaK\'];
    if ~exist(fullfile(outbase))
        mkdir(fullfile(outbase))
    end
    behav_dest = ['Z:\Analysis\motorizedWheel_Analysis\airpuff\behavioral_analysis\' days{iexp} '\'];  
    threshold = -4;
    filename = dir([image_analysis_dest sessions{iexp} '_decon' '*' num2str(threshold) '_TCave_cl.mat']);
    TCave_cl = load([image_analysis_dest filename.name]);
    tc_avg = TCave_cl.TCave_cl;
    filename2 = dir([image_analysis_dest 'getTC\' '*' 'thresh97_coor0.7_mask3D_final.mat']);
    mask3D = load([image_analysis_dest 'getTC\' filename2.name]);
    mask3D = mask3D.mask3D;
   
    %transform TC to dfof
    meanF = mean(tc_avg, 1); %average fluorescence across time for each cell
    df_f = bsxfun(@minus, tc_avg, meanF); % raw F - mean F
    df_f = bsxfun(@rdivide, df_f, meanF); % normalize by meanF
    max_dff = max(df_f, [], 1); %find brightest value of each cell
    df_f_norm = bsxfun(@rdivide, df_f, max_dff); %normalized by brightest value
    
    mask_cell = zeros(size(mask3D(:,:,1)));
    for i = 1:size(tc_avg,2)
        mask_cell(find(mask3D(:,:,i))) = i; %pixels belong to a dendrite is assigned as dendrite n 
    end
    
%     behav_output = load([behav_dest days{iexp} '_behavAnalysis.mat']);
%     frm_stay_cell = behav_output.frames_stay_cell;
%     frm_stay = cell2mat(frm_stay_cell);
%     airpuffon = behav_output.airpuffon1;
%     % find the frames during stationary without airpuff (get rid of frames 300ms after airpuff onset)
%     % 1.stationary without airpuff
%     airpuff_period = [];
%     for a = 1:length(airpuffon)
%         airpuff_period = cat(2,airpuff_period,airpuffon(a):airpuffon(a)+9);%300ms after airpuff onset is 10ish frames
%     end
%     stay_noairpuff = setdiff(frm_stay,airpuff_period);% find the frames in frame_stay but not in airpuff_period
%     
%     df_f_norm_spont = df_f_norm(stay_noairpuff,:);
    
    % no seeds
    sc = 0;
%     % code for making seeds
%     h = figure;
%     imagesc(mask_cell)
%     [i j] = getpts(h);% # of points = # of regions you think
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
            r_dend(id,id2) = triu2vec(corrcoef(tc_avg(:,id),tc_avg(:,id2))); %uppper triangular part of the matrix
        end
    end
%     
%     figure;
%     for i = 2:7
%     idx = kmeans(r_dend,17);
%     [x idx_sort] = sort(idx,'ascend');
%     ind = find(diff(x));
%     r_dend_sort = r_dend(idx_sort,idx_sort);
%     subplot(3,2,i-1)
%     imagesc(r_dend_sort)
%     hline(ind+.5)
%     vline(ind+.5)
%     end
%     
%     print(fullfile(outbase, [days{iexp} '_sortMatrix.pdf']),'-dpdf', '-bestfit')
    
    % number of clusters
    n = 4;
    
    % run meta k means
    [allclusters, allclusters_init, centroidcorr, dendmem, dunnsinitial] = meta_k_means_LG(df_f_norm', 'correlation', sc, n);
    
    %plot output of meta k means
    figure;
    mask_cluster = zeros(size(mask_cell));
    [n1 n2] = subplotn(size(allclusters,1)+2);
    for i = 1:size(allclusters,1)
        mask_cluster_temp = zeros(size(mask_cell));
        for ii = 1:size(allclusters{i,1},2)
            mask_cluster(find(mask_cell ==  allclusters{i,1}(ii))) = i;
            mask_cluster_temp(find(mask_cell ==  allclusters{i,1}(ii))) = 1;
        end
        subplot(n1,n2,i);
        imagesc(mask_cluster_temp);
    end
    
    subplot(n1,n2,i+1); imagesc(mask_cluster);
    print(fullfile(outbase, [days{iexp}, '_clusters_n' num2str(n) '_wholeTC']),'-dpdf','-bestfit');
    
    c_order = [];
    b = 0.5;
    for ic = 1:size(allclusters,1)
        c_order = [c_order allclusters{ic,1}];
        b = [b b(end)+length(allclusters{ic,1})];
    end
    figure; imagesc(r_dend(c_order,c_order))
    hline(b);
    vline(b);
    print(fullfile(outbase, [days{iexp} '_corrclusters_n' num2str(n) '_wholeTC']),'-dpdf','-bestfit');
    
    save(fullfile(outbase, [days{iexp} '_clusters_n' num2str(n) '_wholeTC.mat']), 'allclusters', ...
        'allclusters_init','centroidcorr', 'dendmem', 'dunnsinitial', 'mask_cell', 'mask_cluster', 'sc');
    
    %{
    % if using seeds
    if size(scs,2)>1
        subplot(n,n2,i+2); imagesc(sc_mask)
    end
    print(fullfile(outbase, days{iexp}, [days{iexp} '_clusters_sc.pdf']),'-dpdf')
    %}
    
    
    %plot airpuff response of different clusters
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvFbtm_airpuff_stay_cells = dfOvF_strct.dfOvFbtm_airpuff_stay_cells;
    cluster1_resp = dfOvFbtm_airpuff_stay_cells(:,allclusters{2,1});
    cluster2_resp = dfOvFbtm_airpuff_stay_cells(:,allclusters{3,1});
    x = (1:size(cluster1_resp,1))/30;
    figure;
    subplot(2,1,1);
    fast_errbar(x,cluster1_resp,2,'shaded', true,'color',[0.9686 0.5059 0.7490]);
    xlim([0 1]); ylim([0 0.7]);
    set(gca,'Fontsize',18);
    ylabel('df/F');
    vline(16/30,'k');
    title('cluster 2');
    subplot(2,1,2);
    fast_errbar(x,cluster2_resp,2,'shaded', true,'color',[0.9686 0.5059 0.7490]);
    xlim([0 1]); ylim([0 0.7]);ylabel('df/F');xlabel('time(s)');
    set(gca,'Fontsize',18);
    vline(16/30,'k');
    title('cluster 3');
        
end
