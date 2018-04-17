function [h, p, stats]=scatter_plot(mouseID, y1, y2, col_mat)

cid = 1;
for id = 1:size(mouseID,2)
    %     y1_size = size(y1{id}, 2);
    %     y2_size = size(y2{id}, 2);
    y1_temp = y1{id};
    y1_size = max(size(y1_temp(~isnan(y1_temp))));
    y2_temp = y2{id};
    y2_size = max(size(y2_temp(~isnan(y2_temp))));
    y1_mean = nanmean(y1{id});
    y2_mean = nanmean(y2{id});
    y1_std = nanstd(y1{id});
    y2_std = nanstd(y2{id});
    %     y1_mean = mean(y1{id}, 2);
    %     y2_mean = mean(y2{id},2);
    %     y1_std = std(y1{id},[],2);
    %     y2_std = std(y2{id},[],2);
    hold on
    
    if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
        
        cid = cid + 1;
        
    end
    
    if size(col_mat,1) > 1 && mean(col_mat(1,:) ~= col_mat(2,:))
        id
        if y1_size ~= 1
            errorbarxy_2p(y1_mean, y2_mean,  y1_std./sqrt(y1_size), y2_std./sqrt(y2_size),{'o', col_mat(cid,:), col_mat(cid,:),col_mat(cid,:)})
            hold on
%             scatter(y1_mean, y2_mean, 'MarkerFaceColor', col_mat(cid,:), 'MarkerEdgeColor',col_mat(cid,:), 'LineWidth', 1)
            scatter(y1_mean, y2_mean, 'MarkerEdgeColor',col_mat(cid,:), 'LineWidth', 1)
        else
            scatter(y1_mean, y2_mean, 25, col_mat(cid,:), 'filled');
            hold on
        end
    elseif size(col_mat,1) == 1
%         errorbarxy(y1_mean, y2_mean,  y1_std./sqrt(y1_size), y2_std./sqrt(y2_size),{'o', col_mat(1,:), col_mat(1,:),col_mat(1,:)})
        hold on
        scatter(y1{id}, y2{id}, 4, 'MarkerEdgeColor',col_mat(1,:))
    else
         errorbarxy(y1_mean, y2_mean,  y1_std./sqrt(y1_size), y2_std./sqrt(y2_size),{'o', col_mat(cid,:), col_mat(cid,:),col_mat(cid,:)})
            hold on
            scatter(y1_mean, y2_mean, 'MarkerFaceColor', col_mat(cid,:), 'MarkerEdgeColor',col_mat(cid,:), 'LineWidth', 1)
    end
    
    
end
if size(col_mat,1) == 1%size(col_mat,1) > 1 && mean(col_mat(1,:) == col_mat(2,:)) %size(col_mat,1) == 1 %
    y1 = cell2mat(y1);
    y1_size = max(size(y1));
    y1_mean = nanmean(y1);
    y2 = cell2mat(y2);
    y2_size = max(size(y2));
    y2_mean = nanmean(y2);
    col = [0.9 0 0];
    if y1_std == 0
        y2_std = std(y2);
    end
    h = errorbarxy(y1_mean, y2_mean,  y1_std./sqrt(y1_size), y2_std./sqrt(y2_size),{'o', col, col, col});
    set(h.hMain,'LineWidth', 1);
    [h,p,ci,stats] = ttest(y1, y2);
    
end

end
