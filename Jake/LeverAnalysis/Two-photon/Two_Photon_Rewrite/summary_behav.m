

file_info
out_base = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis\');
for id = 1:size(mouseID,2)
    for rID  = 1:2
        dest_sub  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
        if exist(dest_sub)
%             winopen([dest_sub '_cell_resp_amplitude.pdf']);
            load([dest_sub '_release_movies.mat'])
            load([dest_sub '_press_movies.mat'])
            load([dest_sub '_cell_categories.mat'])
%             load([dest_sub '_cell_resp_amp.mat'])
            load([dest_sub '_cell_resp.mat']);
            stats_s = size(success_movie,1);
            stats_f = size(fail_movie,1);
            RS_cells{id} = unique([release_resp_cells success_resp_cells fail_resp_cells press_resp_cells]);
            stats(id) = stats_s/(stats_s + stats_f);
            seg_cells(id) = size(success_movie,2);
            release_resp_cells_tot(id) = length(release_resp_cells);
            tot_resp(id) = length(unique([success_resp_cells fail_resp_cells press_resp_cells]));
            release_resp_cell = unique([success_resp_cells fail_resp_cells]);
            tot_release(id) = length(release_resp_cell);
            pct_release(id) = tot_release(id)./tot_resp(id);
%             release_resp_amp(id) = mean(mean(release_peak));
            release_resp_RS = mean((release_resp(:,RS_cells{id})-release_base(:,RS_cells{id})),1);
            release_resp_RS_mean(id) = mean(release_resp_RS,2);
        end
    end
end

col_mat = [ 0.9  0.9  0;
    1  0  1;
    0  1  1;
    0.5  0  0;
    0  1  0;
    0  0  1;
    1  0.6  1;
    0  0  0;
    1  0.8 0.4
    0  0.5 0.7
    0.5 0.4 0];

fig=figure;
subplot(2,2,1)
cid = 1;
plotMouse = unique(mouseID);
for id = 1:size(mouseID,2)
    
    if id > 1
        if strcmp(mouseID{id},mouseID{id-1})
            scatter(seg_cells(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        else
            cid = cid + 1; 
            scatter(seg_cells(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
            
            
        end
    else
        scatter(seg_cells(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        
    end
    hold on
end

x = 1:0.79:80;
y = 0:0.01:1;
xlim([0 150])
ylim([0 1])
% plot(x,y,'-k')
ylabel('correct %')
xlabel('# of segmented cells')

subplot(2,2,2)
cid = 1;
for id = 1:size(mouseID,2)
    
    if id > 1
        if strcmp(mouseID{id},mouseID{id-1})
            scatter(release_resp_cells_tot(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        else
            cid = cid + 1; 
            scatter(release_resp_cells_tot(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
            
        end
    else
        scatter(release_resp_cells_tot(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
    end
    hold on
end
x=1:0.79:80;
y = 0:0.01:1;
% plot(x,y,'-k')
ylabel('correct %')
xlabel('# of release resp cells')

subplot(2,2,3)
cid = 1;
for id = 1:size(mouseID,2)
    
    if id > 1
        if strcmp(mouseID{id},mouseID{id-1})
            scatter(pct_release(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        else
            cid = cid + 1; 
            scatter(pct_release(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
            
        end
    else
        scatter(pct_release(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
    end
    hold on
end
x = 0:.1:1;
y = x;
xlim([0 1.1]);
% plot(x,y,'-k')
ylabel('correct %')
xlabel('fraction of release resp cells')

subplot(2,2,4)
cid = 1;
for id = 1:size(mouseID,2)
    
    if id > 1
        if strcmp(mouseID{id},mouseID{id-1})
            scatter(release_resp_RS_mean(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
        else
            cid = cid + 1; 
            scatter(release_resp_RS_mean(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
            
        end
    else
        scatter(release_resp_RS_mean(id), stats(id), 'MarkerEdgeColor',col_mat(cid,:),'MarkerFaceColor',col_mat(cid,:))
    end
    hold on
end
x=0:0.004:0.4;
y = 0:0.01:1;
xlim([0 0.4])
% plot(x,y,'-k')
ylabel('correct %')
xlabel('release resp amp')
legend(mouseID)
saveas(fig, [out_base, 'summary_behav.fig']);
print([out_base 'Summary_behav_stats.eps'], '-depsc');
print([out_base 'Summary_behav_stats.pdf'], '-dpdf');