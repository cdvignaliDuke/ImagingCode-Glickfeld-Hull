clear;
%NEED TO UPDATE THIS SO IT ACCESSES SPREADSHEET INSTEAD OF JUST WRITING IN THE NAMES
sessions = {'191114_img1040','191115_img1039','191115_img1041','191115_img1042','200316_img1064_airpuff_2'};
image_analysis_base  = 'Z:\Analysis\Airpuff_analysis\imaging_analysis\';%stores the data on crash in the movingDots analysis folder
color_code = {'c','r','y','g'};

%% resp amp stay vs. running
for ii = 4
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvFbtm_airpuff_stay_cells = dfOvF_strct.dfOvFbtm_airpuff_stay_cells;
    dfOvFbtm_airpuff_run_cells = dfOvF_strct.dfOvFbtm_airpuff_allrun_cells;
    baseline_cell_stay = zeros(1,size(dfOvFbtm_airpuff_stay_cells,2));
    resp_cell_stay = zeros(1,size(dfOvFbtm_airpuff_stay_cells,2));
    for c = 1: size(dfOvFbtm_airpuff_stay_cells,2)
        % calculate baseline F for each cell
        temp = dfOvFbtm_airpuff_stay_cells(:,c);
        baseline_cell_stay(c) = mean(temp(1:15));
        resp_cell_stay(c) = max(temp(17:22)) - baseline_cell_stay(c);
    end
    
    baseline_cell_run = zeros(1,size(dfOvFbtm_airpuff_stay_cells,2));
    resp_cell_run = zeros(1,size(dfOvFbtm_airpuff_stay_cells,2));
    for c = 1: size(dfOvFbtm_airpuff_run_cells,2)
        % calculate baseline F for each cell
        temp = dfOvFbtm_airpuff_run_cells(:,c);
        baseline_cell_run(c) = mean(temp(1:15));
        resp_cell_run(c) = max(temp(17:22)) - baseline_cell_run(c);
    end
    
    min_resp_run = min(resp_cell_run);
    min_resp_stay = min(resp_cell_stay);
    max_resp_run = max(resp_cell_run);
    max_resp_stay = max(resp_cell_stay);
    
    filename2 = dir([image_analysis_dest 'getTC\' '*' 'thresh97.5_coor0.8_mask3D_final.mat']);
    mask3D = load([image_analysis_dest 'getTC\' filename2.name]);
    mask3D = mask3D.mask3D;
    
    mask_cell_stay = zeros(size(mask3D(:,:,1)));
    mask_cell_run = zeros(size(mask3D(:,:,1)));
    for c = 1:length(resp_cell_stay)
        mask_cell_stay = mask_cell_stay + mask3D(:,:,c).*resp_cell_stay(c);
        mask_cell_run = mask_cell_run + mask3D(:,:,c).*resp_cell_run(c);
    end
    isc1 = figure; 
    imagesc(mask_cell_stay,[min(min_resp_run,min_resp_stay) max(max_resp_run,max_resp_stay)]); colorbar;
    title('stationary response amplitude');
    set(gca,'Fontsize',8);
    isc1.Units = 'centimeters';
    isc1.Position = [1 3 12 5];
    fig_name = ['airpuff_resp_stayamp_mask' sessions{ii}];
    path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
    savefig([path,fig_name]);
    %savefig([image_analysis_dest 'resp\' sessions{ii} '_resp_amp_mask']);

    
    isc2 = figure;
    imagesc(mask_cell_run,[min(min_resp_run,min_resp_stay) max(max_resp_run,max_resp_stay)]);colorbar;% make the scale of the 2 behavioral states the same: the smallest and biggest value among resp_prob
    title('running response amplitude');
    set(gca,'Fontsize',8);
    isc2.Units = 'centimeters';
    isc2.Position = [1 3 12 5];
    fig_name = ['airpuff_resp_runamp_mask' sessions{ii}];
    path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
    savefig([path,fig_name]);
%     print(isc2,[path,fig_name],'-r600','-dpdf');
%    savefig([image_analysis_dest 'resp\' sessions{ii} '_resp_amp_mask']);
    
    
%     figure;
%     scatter(resp_cell_stay,resp_cell_run,'filled','MarkerEdgeColor',...
%         [0.6000 0.6000 0.6000],'MarkerFaceColor',[0.6000 0.6000 0.6000]);
%     axis square; xlabel('stay');ylabel('run');
%     line = refline(1,0);
%     line.Color = [0.9686 0.5059 0.7490];
%     line.LineWidth = 1.5;
%     set(gca,'Fontsize',18);
%     xlim([-0.2 1.2]);ylim([-0.2 1.2]);
%     title(['airpuff response amplitude ' sessions{ii}]);
%     savefig([image_analysis_dest '\' sessions{ii} '_airpuff_resp_amp_scatter']);
%     save([image_analysis_dest '\' sessions{ii} '_dfOvF.mat' ],...
%         'baseline_cell_stay','resp_cell_stay','baseline_cell_run',...
%         'resp_cell_run','-append');
end

%% modulation index (running resp - stay resp/stay resp)
for ii = 1: length(sessions)
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    %plot all of the trials of all cells into 1 histogram
    
    resp_cell_stay = dfOvF_strct.resp_cell_stay;
    resp_cell_run = dfOvF_strct.resp_cell_run;
    resp_diff = resp_cell_run - resp_cell_stay;
    modu_inx = resp_diff./abs(resp_cell_stay);
    
    filename2 = dir([image_analysis_dest 'getTC\' '*' 'thresh97.5_coor0.8_mask3D.mat']);
    mask3D = load([image_analysis_dest 'getTC\' filename2.name]);
    mask3D = mask3D.mask3D;
    
    mask_cell = zeros(size(mask3D(:,:,1)));
    for c = 1:length(resp_cell_run)
        mask_cell = mask_cell + mask3D(:,:,c).*modu_inx(c);
    end
    figure; imagesc(mask_cell,[min(modu_inx) 1.5]);
    %clims([min(modu_inx) 1]);
    title([sessions{ii} 'airpuff modulation index of individual cells']);
    savefig([image_analysis_dest 'resp\' '_respMI_mask']);
    
    save([image_analysis_dest '\' sessions{ii} '_dfOvF.mat' ],...
        'modu_inx','-append');
end


%% resp prob
for ii = 1:4
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvFbtm_airpuff_stay_mat = dfOvF_strct.dfOvFbtm_airpuff_stay_mat;      % trial*frame*cell
    dfOvFbtm_airpuff_allrun_mat = dfOvF_strct.dfOvFbtm_airpuff_allrun_mat;
    % test if the airpuff response is significant
    h_stay = zeros(1,size(dfOvFbtm_airpuff_stay_mat,3));
    p_stay = zeros(1,size(dfOvFbtm_airpuff_stay_mat,3));
    base_trials_stay = zeros(size(dfOvFbtm_airpuff_stay_mat,1),size(dfOvFbtm_airpuff_stay_mat,3));% trails*cells
    stdbase_trials_stay = zeros(size(dfOvFbtm_airpuff_stay_mat,1),size(dfOvFbtm_airpuff_stay_mat,3));
    max_trials_stay = zeros(size(dfOvFbtm_airpuff_stay_mat,1),size(dfOvFbtm_airpuff_stay_mat,3));
    resp_trials_stay = zeros(size(dfOvFbtm_airpuff_stay_mat,1),size(dfOvFbtm_airpuff_stay_mat,3));%trials*cells
    resp_period_ave_stay = zeros(size(dfOvFbtm_airpuff_stay_mat,1),size(dfOvFbtm_airpuff_stay_mat,3));%trials*cells
    h_run = zeros(1,size(dfOvFbtm_airpuff_allrun_mat,3));
    p_run = zeros(1,size(dfOvFbtm_airpuff_allrun_mat,3));
    base_trials_run = zeros(size(dfOvFbtm_airpuff_allrun_mat,1),size(dfOvFbtm_airpuff_allrun_mat,3));
    stdbase_trials_run = zeros(size(dfOvFbtm_airpuff_allrun_mat,1),size(dfOvFbtm_airpuff_allrun_mat,3));
    max_trials_run = zeros(size(dfOvFbtm_airpuff_allrun_mat,1),size(dfOvFbtm_airpuff_allrun_mat,3));
    resp_trials_run = zeros(size(dfOvFbtm_airpuff_allrun_mat,1),size(dfOvFbtm_airpuff_allrun_mat,3));
    resp_period_ave_run = zeros(size(dfOvFbtm_airpuff_allrun_mat,1),size(dfOvFbtm_airpuff_allrun_mat,3));
%     [n1,n2] = subplotn(size(dfOvFbtm_airpuff_allrun_mat,3)+2);
%     [n3,n4] = subplotn(size(dfOvFbtm_airpuff_allrun_mat,3)+2);
    
    for c = 1: size(dfOvFbtm_airpuff_stay_mat,3)
        temp1 = dfOvFbtm_airpuff_stay_mat(:,1:15,c);
        base_trials_stay(:,c) = mean(temp1,2);
        stdbase_trials_stay(:,c) = std(temp1,0,2);
        temp2 = dfOvFbtm_airpuff_stay_mat(:,17:22,c);
        max_trials_stay(:,c) = max(temp2,[],2);
        resp_trials_stay(:,c) = max_trials_stay(:,c) - base_trials_stay(:,c);
        resp_period_ave_stay(:,c) = mean(temp2,2);
        [h_stay(c),p_stay(c)] = ttest(resp_period_ave_stay(:,c),base_trials_stay(:,c),'Tail','right');
        
        temp3 = dfOvFbtm_airpuff_allrun_mat(:,1:15,c);
        base_trials_run(:,c) = mean(temp3,2);
        stdbase_trials_run(:,c) = std(temp3,0,2);
        temp4 = dfOvFbtm_airpuff_allrun_mat(:,17:22,c);
        max_trials_run(:,c) = max(temp4,[],2);
        resp_trials_run(:,c) = max_trials_run(:,c) - base_trials_run(:,c);
        resp_period_ave_run(:,c) = mean(temp4,2);
        [h_run(c),p_run(c)] = ttest(resp_period_ave_run(:,c),base_trials_run(:,c),'Tail','right');
        
        %make histograms of responses during stationary and running for each cell
%         fig_resp_staycells = figure(1);
%         hold on;
%         subplot(n1,n2,c);
%         hist(resp_trials_stay(:,c)); hold on;
%         fig_resp_runcells = figure(2);
%         hold on;
%         subplot(n3,n4,c);
%         hist(resp_trials_run(:,c)); hold on;
    end
    %     figure(1); supertitle(['resp_amp_stay' sessions{ii}]);
    %     figure(2); supertitle(['resp_amp_run' sessions{ii}]);
    %     saveas(fig_resp_staycells,[image_analysis_dest '\' sessions{ii} '_respAmp_cells_stay'])
    %     saveas(fig_resp_runcells, [image_analysis_dest '\' sessions{ii} '_respAmp_cells_run']);
    
    %reshape the trial*cell matrix into vector and plot histogram. the histogram for each cell is hard to visualize when determining response threshold
%     resp_trials_run_vec = reshape(resp_trials_run,[],1);
%     resp_trials_stay_vec = reshape(resp_trials_stay,[],1);
%     stdbase_trials_stay_vec = reshape(stdbase_trials_stay,[],1);
%     stdbase_trials_run_vec = reshape(stdbase_trials_run,[],1);
%     resp_hist = figure;
%     subplot(2,2,1);histogram(resp_trials_run_vec,'BinWidth',0.05);title('resp amp run');
%     subplot(2,2,2);histogram(resp_trials_stay_vec,'BinWidth',0.05);title('resp amp stay');
%     subplot(2,2,3);histogram(stdbase_trials_run_vec,'BinWidth',0.05);title('std run');
%     subplot(2,2,4);histogram(stdbase_trials_stay_vec,'BinWidth',0.05);title('std stay');
%     supertitle(sessions{ii});
%     saveas(resp_hist,[image_analysis_dest 'resp\' sessions{ii} '_resp_and_std_hist']);
    
%     %based on the histogram, choose a response threshold and calulate the
%     %response probablity for each cell
%     resp_thres = 0.2; %an increase of 0.2 in df/F
    ntrials_stay = size(resp_trials_stay,1);
    ntrials_run = size(resp_trials_run,1);
%     resp_prob_run = zeros(1,size(dfOvFbtm_airpuff_stay_mat,3));
%     resp_prob_stay = zeros(1,size(dfOvFbtm_airpuff_stay_mat,3));
    resp_prob_run_2std = zeros(1,size(dfOvFbtm_airpuff_stay_mat,3));
    resp_prob_stay_2std = zeros(1,size(dfOvFbtm_airpuff_stay_mat,3));
    resp_prob_run_3std = zeros(1,size(dfOvFbtm_airpuff_stay_mat,3));
    resp_prob_stay_3std = zeros(1,size(dfOvFbtm_airpuff_stay_mat,3));
    for c = 1: size(dfOvFbtm_airpuff_stay_mat,3)
        % use ave baseline+2std
        resp_prob_run_2std(c) = sum(max_trials_run(:,c)> base_trials_run(:,c)+2*stdbase_trials_run(:,c))/ntrials_run; %when you get very few numbers of running trials, the response probability during running of a bunch of cells will be the same
        resp_prob_stay_2std(c) = sum(max_trials_stay(:,c)>base_trials_stay(:,c)+2*stdbase_trials_stay(:,c))/ntrials_stay;
        resp_prob_run_3std(c) = sum(max_trials_run(:,c)> base_trials_run(:,c)+3*stdbase_trials_run(:,c))/ntrials_run; %when you get very few numbers of running trials, the response probability during running of a bunch of cells will be the same
        resp_prob_stay_3std(c) = sum(max_trials_stay(:,c)>base_trials_stay(:,c)+3*stdbase_trials_stay(:,c))/ntrials_stay;
%         % use a solid threshold
%         resp_prob_run(c) = sum(resp_trials_run(:,c)>resp_thres)/ntrials_run;
%         resp_prob_stay(c) = sum(resp_trials_stay(:,c)>resp_thres)/ntrials_stay;
    end
    
%     mean_resp_prob_run_2std = mean(resp_prob_run_2std);
%     mean_resp_prob_stay_2std = mean(resp_prob_stay_2std);
    
%     scatter_prob = figure;
%     scatter(resp_prob_stay_2std,resp_prob_run_2std,8,'filled','MarkerEdgeColor',...
%         [0.45 0.45 0.45],'MarkerFaceColor',[0.45 0.45 0.45]); hold on;
%     scatter(mean_resp_prob_stay_2std,mean_resp_prob_run_2std,8,'filled','MarkerEdgeColor',...
%         [0.8431 0.0980 0.1098],'MarkerFaceColor',[0.8431 0.0980 0.1098]); hold on;
%     axis square; xlabel('response probability (stationary)');ylabel('response probability (running)');
%     xlim([0 1]);ylim([0 1]);
%     line = refline(1,0);
%     line.Color = 'r';
%     line.LineWidth = 1;
%     set(gca,'Fontsize',8);
%     scatter_prob.Units = 'centimeters';
%     scatter_prob.Position = [1 3 5 5];
%     %fig_name = ['airpuff_resp_prob_scatter_' sessions{ii}];
%     %path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
%     %print(scatter_prob,[path,fig_name],'-r600','-depsc');
%     title(['airpuff response probability, 3std' sessions{ii}]);
%     savefig([image_analysis_dest 'resp\' sessions{ii} '_airpuff_resp_prob_3std_scatter']);
    
%     filename2 = dir([image_analysis_dest 'getTC\' '*' 'thresh97_coor0.8_mask3D.mat']);
%     mask3D = load([image_analysis_dest 'getTC\' filename2.name]);
%     mask3D = mask3D.mask3D;
%     
%     mask_cell_stay = zeros(size(mask3D(:,:,1)));
%     mask_cell_run = zeros(size(mask3D(:,:,1)));
%     for c = 1:length(resp_prob_stay_std)
%         mask_cell_stay = mask_cell_stay + mask3D(:,:,c).*resp_prob_stay_std(c);
%         mask_cell_run = mask_cell_run + mask3D(:,:,c).*resp_prob_run_std(c);
%     end
%     isc1 = figure; 
%     imagesc(mask_cell_stay,[0,1]); colorbar;
%     title('stationary');
%     set(gca,'Fontsize',8);
%     isc1.Units = 'centimeters';
%     isc1.Position = [1 3 12 5];
%     fig_name = ['airpuff_resp_stayprob_mask' sessions{ii}];
%     path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
%     savefig([path,fig_name]);
%     
%     isc2 = figure;
%     imagesc(mask_cell_run,[0,1]);colorbar;% make the scale of the 2 behavioral states the same: the smallest and biggest value among resp_prob
%     title('running');
%     set(gca,'Fontsize',8);
%     isc2.Units = 'centimeters';
%     isc2.Position = [1 3 12 5];
%     fig_name = ['airpuff_resp_runprob_mask' sessions{ii}];
%     path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
%     savefig([path,fig_name]);
    %print(isc2,[path,fig_name],'-r600','-dpdf');
%     savefig([image_analysis_dest 'resp\' sessions{ii} '_resp_prob_std_mask']);

    save([image_analysis_dest '\' sessions{ii} '_dfOvF.mat' ],...
        'base_trials_stay','base_trials_run',...
        'stdbase_trials_run','stdbase_trials_stay',...
        'max_trials_run','max_trials_stay',...
        'resp_trials_stay','resp_trials_run',...
        'resp_prob_run_2std','resp_prob_stay_2std',...
        'resp_prob_run_3std','resp_prob_stay_3std',...
        'h_run','h_stay','p_run','p_stay','-append');
 %       'resp_prob_run','resp_prob_stay',...
        
end

%% color code response amlitude plot: do the cells that has a hight amp also has a high prob?
% need: response amplitude, and result of response probablity t-test
for ii = 3
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    resp_cell_stay = dfOvF_strct.resp_cell_stay;
    resp_cell_run = dfOvF_strct.resp_cell_run;
    h_stay = dfOvF_strct.h_stay;
    h_run = dfOvF_strct.h_run;
    
    colorcode = {[0.8941    0.1020    0.1098],[0.2157    0.4941    0.7216],[0.3020    0.6863    0.2902],[0.4 0.4 0.4]};
    amp_color = figure;hold on;
    scatter(100,100,'filled','MarkerEdgeColor',colorcode{1},'MarkerFaceColor',colorcode{1});
    scatter(100,100,'filled','MarkerEdgeColor',colorcode{2},'MarkerFaceColor',colorcode{2});
    scatter(100,100,'filled','MarkerEdgeColor',colorcode{3},'MarkerFaceColor',colorcode{3});
    scatter(100,100,'filled','MarkerEdgeColor',colorcode{4},'MarkerFaceColor',colorcode{4});
    for c = 1:length(h_stay)
        if h_stay(c) == 1&& h_run(c) == 0
            scatter(resp_cell_stay(c),resp_cell_run(c),8,'filled','MarkerEdgeColor',...
                colorcode{1},'MarkerFaceColor',colorcode{1});
        elseif h_stay(c) == 0 && h_run(c) == 1
            scatter(resp_cell_stay(c),resp_cell_run(c),8,'filled','MarkerEdgeColor',...
                colorcode{2},'MarkerFaceColor',colorcode{2});
        elseif h_run(c) == 1 && h_stay(c) == 1
            scatter(resp_cell_stay(c),resp_cell_run(c),8,'filled','MarkerEdgeColor',...
                colorcode{3},'MarkerFaceColor',colorcode{3});
        elseif h_run(c) == 0 && h_stay(c) == 0
            scatter(resp_cell_stay(c),resp_cell_run(c),8,'filled','MarkerEdgeColor',...
                colorcode{4},'MarkerFaceColor',colorcode{4});
        end
    end
    axis square; xlabel('stay');ylabel('run');
    xlim([-0.2 1]);ylim([-0.2 1]);
    line = refline(1,0);
    line.Color = [0 0 0];
    line.LineWidth = 1;
    %legend('responsive only stationary','responsive only running','responsive both','responsive neither'); legend('boxoff');
    %title(sessions{ii})
    hold off;
    amp_color.Units = 'centimeters';
    amp_color.Position = [3 3 5 5];
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',8);
    fig_name = ['resp_amp_colored',sessions{ii}];
    %path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
    path = [image_analysis_dest 'resp\'];
    orient(amp_color,'landscape')
    print(amp_color,[path,fig_name],'-r600','-depsc');
    
end

%% response amplitude: pull out ONLY trials that has a response
for ii = 1:4
    image_analysis_dest = [image_analysis_base, sessions{ii}, '\'];
    dfOvF_strct = load([image_analysis_dest sessions{ii} '_dfOvF.mat']);
    dfOvFbtm_airpuff_stay_mat = dfOvF_strct.dfOvFbtm_airpuff_stay_mat;      % trial*frame*cell
    dfOvFbtm_airpuff_allrun_mat = dfOvF_strct.dfOvFbtm_airpuff_allrun_mat;
    base_trials_stay = dfOvF_strct.base_trials_stay;
    base_trials_run = dfOvF_strct.base_trials_run;
    stdbase_trials_run = dfOvF_strct.stdbase_trials_run;
    stdbase_trials_stay = dfOvF_strct.stdbase_trials_stay;
    max_trials_run = dfOvF_strct.max_trials_run;
    max_trials_stay = dfOvF_strct.max_trials_stay;
    h_stay = dfOvF_strct.h_stay;
    h_run = dfOvF_strct.h_run;
    resp_trials_stay = dfOvF_strct.resp_trials_stay; %trials*cells
    resp_trials_run = dfOvF_strct.resp_trials_run;
    
    % pull out responsive trials for each cell
    cell_resp_run = cell(1,size(dfOvFbtm_airpuff_stay_mat,3));
    cell_resp_stay = cell(1,size(dfOvFbtm_airpuff_stay_mat,3));
    for c = 1: size(dfOvFbtm_airpuff_stay_mat,3)
        cell_resp_stay{c} = [];
        cell_resp_run{c} = [];
        for t = 1: size(dfOvFbtm_airpuff_stay_mat,1)
            if max_trials_stay(t,c)> base_trials_stay(t,c)+3*stdbase_trials_stay(t,c)
                cell_resp_stay{c} = cat(1,cell_resp_stay{c},resp_trials_stay(t,c));
            end
        end
        for t = 1: size(dfOvFbtm_airpuff_allrun_mat,1)
            if max_trials_run(t,c)> base_trials_run(t,c)+3*stdbase_trials_run(t,c)
                cell_resp_run{c} = cat(1,cell_resp_run{c},resp_trials_run(t,c));
            end
        end
    end
    % average across responsive trials of each cell
    resp_stay_ONLYh1 = cellfun(@mean,cell_resp_stay);
    resp_run_ONLYh1 = cellfun(@mean,cell_resp_run);
    % if some cells don't have any responsive trials, make the Nan into zero
    resp_stay_ONLYh1(isnan(resp_stay_ONLYh1)) = 0;
    resp_run_ONLYh1(isnan(resp_run_ONLYh1)) = 0;
    
%     colorcode = {[0.8941    0.1020    0.1098],[0.2157    0.4941    0.7216],[0.3020    0.6863    0.2902],[0.4 0.4 0.4]};
%     amp_color_onlyh1 = figure;hold on;
%     scatter(100,100,'filled','MarkerEdgeColor',colorcode{1},'MarkerFaceColor',colorcode{1});
%     scatter(100,100,'filled','MarkerEdgeColor',colorcode{2},'MarkerFaceColor',colorcode{2});
%     scatter(100,100,'filled','MarkerEdgeColor',colorcode{3},'MarkerFaceColor',colorcode{3});
%     scatter(100,100,'filled','MarkerEdgeColor',colorcode{4},'MarkerFaceColor',colorcode{4});
%     for c = 1:length(h_stay)
%         if h_stay(c) == 1&& h_run(c) == 0
%             scatter(resp_stay_ONLYh1(c),resp_run_ONLYh1(c),8,'filled','MarkerEdgeColor',...
%                 colorcode{1},'MarkerFaceColor',colorcode{1});
%         elseif h_stay(c) == 0 && h_run(c) == 1
%             scatter(resp_stay_ONLYh1(c),resp_run_ONLYh1(c),8,'filled','MarkerEdgeColor',...
%                 colorcode{2},'MarkerFaceColor',colorcode{2});
%         elseif h_run(c) == 1 && h_stay(c) == 1
%             scatter(resp_stay_ONLYh1(c),resp_run_ONLYh1(c),8,'filled','MarkerEdgeColor',...
%                 colorcode{3},'MarkerFaceColor',colorcode{3});
%         elseif h_run(c) == 0 && h_stay(c) == 0
%             scatter(resp_stay_ONLYh1(c),resp_run_ONLYh1(c),8,'filled','MarkerEdgeColor',...
%                 colorcode{4},'MarkerFaceColor',colorcode{4});
%         end
%     end
%     % xlim([-0.2 1.2]);ylim([-0.2 1.2]);
%     line = refline(1,0);
%     line.Color = [0 0 0];
%     line.LineWidth = 1;
%     %legend('responsive only stationary','responsive only running','responsive both','responsive neither'); legend('boxoff');
%     axis square; xlabel('stay');ylabel('run');
%     %title(['response amplitude' sessions{ii}]);
%     hold off;
%     
%     amp_color_onlyh1.Units = 'centimeters';
%     amp_color_onlyh1.Position = [3 3 5 5];
%     a = get(gca,'XTickLabel');
%     set(gca,'XTickLabel',a,'FontSize',8);
%     fig_name = ['resp_amp_colored_only_responsive_trials',sessions{ii}];
%     %path = [image_analysis_dest 'resp\'];
%     path = 'Z:\Analysis\figures\figure5_airpuff_resp\';
%     orient(amp_color_onlyh1,'landscape')
%     print(amp_color_onlyh1,[path,fig_name],'-r600','-depsc');
    
     save([image_analysis_dest '\' sessions{ii} '_dfOvF.mat' ],...
        'resp_stay_ONLYh1','resp_run_ONLYh1','-append');

%     dfOvFbtm_airpuff_respstay_mat = []; %frames * trials, these trials will be trials from different cells
%     dfOvFbtm_airpuff_resprun_mat = [];
%     for c = 1: size(dfOvFbtm_airpuff_stay_mat,3) %for each cell
%         for t1 = 1: size(dfOvFbtm_airpuff_stay_mat,1) % for each trial
%             if max_trials_stay(t1,c)> base_trials_stay(t1,c)+2*stdbase_trials_stay(t1,c) % if this trial is responsive
%                dfOvFbtm_airpuff_respstay_mat = cat(1,dfOvFbtm_airpuff_respstay_mat,dfOvFbtm_airpuff_stay_mat(t1,:,c));
%             end
%         end
%         for t1 = 1: size(dfOvFbtm_airpuff_run_mat,1) % for each trial
%             if max_trials_run(t1,c)> base_trials_run(t1,c)+2*stdbase_trials_run(t1,c) % if this trial is responsive
%                dfOvFbtm_airpuff_resprun_mat = cat(1,dfOvFbtm_airpuff_resprun_mat,dfOvFbtm_airpuff_run_mat(t1,:,c));
%             end
%         end
%         
%     end
    
end




