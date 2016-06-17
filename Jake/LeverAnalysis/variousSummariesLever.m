%spatial correlation scatter plot 

img24 = [0.60, 0.81, 0.47];
img25 = [0.93, 0.96, 0.80, 0.92, 0.79, 0.80];
img27 = [0.76, 0.65, 0.64];
%img28 = 

valuesY = [img24, img25, img27];
categoriesX = [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3];

corr = scatter( categoriesX, valuesY, 20, 'b', 'filled')
xlim([0 4]);
ylim([0 1]);
corr.XTick = [0 1 2 3 4];

%% DISTANCE FINDER
for i = 1: cluster25.num_cluster;
xx = reshape(cluster25.roi_mask(i,:),50,50);
[row, column] = find(xx);
avgR = mean(row)
avgC = mean(column)
hold on; plot(avgC,avgR, 'w')
R(i) = avgR;
C(i) = avgC; 
end

d = [];
d = sqrt((C(4)-C(3))^2+(R(4)-R(3))^2)

D_C_img24 = [11.2, 0.61; 9.2, 0.81; 18.9, 0.47]   %3
D_C_img25 = [7.5, 0.93; 9.8, 0.97; 17.2, 0.80; 13.0, 0.92; 12.1, 0.79; 14.8, 0.80]    %6
D_C_img27 = [ 13.0, 0.76; 15.7, 0.65; 10.0, 0.64]    %3

D_C_img24(:,1) = D_C_img24(:,1)/0.01325
D_C_img25(:,1) = D_C_img25(:,1)/0.01325
D_C_img27(:,1) = D_C_img27(:,1)/0.01325

figure;
hold on
scatter(D_C_img24(:,1)', D_C_img24(:,2)', 'or', 'fill')
scatter(D_C_img25(:,1)', D_C_img25(:,2)', 'b', 'fill')
scatter(D_C_img27(:,1)', D_C_img27(:,2)', 'm', 'fill')
xlabel('distance')
ylabel('correlation')
ylim([0 1])
xlim([0 1500])

x = -.1:.1:1;
y = x;
 plot(x,y,'k')
 
 
 
 %% ROI EXAMPLE TRACE

 hold on
 for i = 1:cluster.num_cluster;
    plot(cluster.roi_position{i}(:,1), cluster.roi_position{i}(:,2), 'w',...
        'LineWidth',2)
 end 
 
 
 
 
 
 
 
 
