
color_cell = {'k', 'b', 'g', 'm', 'c', 'r', 'y', 'k', 'b', 'g', 'm', 'c', 'r', 'y', 'k', 'b', 'g', 'm', 'c', 'r', 'y', 'k', 'b', 'g', 'm', 'c', 'r', 'y', 'k', 'b', 'g', 'm', 'c', 'r', 'y', 'k', 'b', 'g', 'm', 'c', 'r', 'y'};
x_axis = ([1:size(NR_movie_nolick, 3)]-61)*33;
NR_movie_nolick = squeeze(nanmean(NR_movie_nolick,1));
OR_movie_nolick = squeeze(nanmean(OR_movie_nolick,1));

for ii = [39,50,53,60];    % 1: size(NR_movie_nolick,1); % [1,2,3,4,11,14,15,19,20,22,23,29,35,36,37,40,41,42,43,50,61];img90  %   % [7,20,30,35,59,71]img91
    figure; hold on;
    plot(x_axis, NR_movie_nolick(ii,:), 'g');
    plot(x_axis, OR_movie_nolick(ii,:), 'r');
    %ylim([-0.04 0.16]); xlim([-1000 2000]);
    vline(0, 'k'); vline(600, 'r');
    title(['img 94 day 1 mask number = ', num2str(ii)]);
    pause 
end

%figure; imagesc(squeeze(mask3D(:,:,20)));