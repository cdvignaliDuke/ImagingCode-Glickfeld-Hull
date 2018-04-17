function corrMatrix = corrLick(success_movie, success_lick)
for ic = 1:size(success_movie,2)
    success_movie_cell = squeeze(success_movie(:,ic,16:end));
   
    if ic == 1
        success_lick = success_lick(:, 16:end);
        check_lick = sum(success_lick, 2);
        success_lick(check_lick == 0, :) = [];
    end
    success_movie_cell(check_lick == 0, :) = [];
    for ie = 1:size(success_movie_cell,1)
        
        [r(ie), p(ie)] = corr(success_movie_cell(ie,:)', success_lick(ie,:)', 'type','Spearman');
    end
    corrMatrix(ic,:) = [mean(r) mean(p)];
end
end