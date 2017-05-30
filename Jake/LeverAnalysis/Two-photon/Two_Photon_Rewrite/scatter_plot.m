function [h, p]=scatter_plot(mouseID, y1, y2, col_mat)

cid = 1;
for id = 1:size(mouseID,2)
    y1_size = size(y1{id}, 2);
    y2_size = size(y2{id}, 2);
    y1_mean = mean(y1{id}, 2);
    y2_mean = mean(y2{id},2);
    y1_std = std(y1{id},[],2);
    y2_std = std(y2{id},[],2);
    hold on
    
        if id > 1 && ~strcmp(mouseID{id},mouseID{id-1})
            
            cid = cid + 1;
             
        end
        
        if col_mat(1,:) ~= col_mat(2,:)
            errorbarxy(y1_mean, y2_mean,  y1_std./sqrt(y1_size), y2_std./sqrt(y2_size),{'o', col_mat(cid,:), col_mat(cid,:),col_mat(cid,:)})
            hold on 
            scatter(y1_mean, y2_mean, 'MarkerEdgeColor',col_mat(cid,:), 'LineWidth', 1)
        else
            scatter(y1{id}, y2{id}, 4, 'MarkerEdgeColor',col_mat(cid,:))
        end
    
    
end
if col_mat(1,:) == col_mat(2,:)
    y1 = cell2mat(y1);
    y1_size = size(y1,2);
    y1_mean = mean(y1,2);
    y2 = cell2mat(y2);
    y2_size = size(y2,2);
    y2_mean = mean(y2,2);
    col = [0.9 0 0];
    h = errorbarxy(y1_mean, y2_mean,  y1_std./sqrt(y1_size), y2_std./sqrt(y2_size),{'o', col, col, col});
    set(h.hMain,'LineWidth', 1);
    [h, p] = ttest(y1, y2);

end

end
