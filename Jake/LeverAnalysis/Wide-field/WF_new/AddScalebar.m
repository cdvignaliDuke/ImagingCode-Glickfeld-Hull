load('cluster_wholeMovie_clean.mat')

% 151009_img30
figure;imagesc(clusters)
hold on; subplot(2,3,1);line([5 81], [250 250], 'color', 'white', 'LineWidth', 2)
subplot(2,3,1);text('position', [230 81], 'fontsize', 10, 'string', '100 µm', 'color', 'white')

% 160314_img38
figure;imagesc(clusters)
hold on; line([5 16], [75 75], 'color', 'white', 'LineWidth', 2)
text('position', [4 71], 'fontsize', 14, 'string', '200 µm', 'color', 'white')

% 160606_img46
figure;imagesc(clusters)
hold on; line([5 16], [60 60], 'color', 'white', 'LineWidth', 2)
text('position', [4 57], 'fontsize', 14, 'string', '200 µm', 'color', 'white')

% 160722_img53
figure;imagesc(clusters)
hold on; line([5 16], [50 50], 'color', 'white', 'LineWidth', 2)
text('position', [5.5 47], 'fontsize', 14, 'string', '200 µm', 'color', 'white')

% 160904_img55
figure;imagesc(clusters)
hold on; line([5 16], [105 105], 'color', 'white', 'LineWidth', 2)
text('position', [4 100], 'fontsize', 14, 'string', '200 µm', 'color', 'white')

