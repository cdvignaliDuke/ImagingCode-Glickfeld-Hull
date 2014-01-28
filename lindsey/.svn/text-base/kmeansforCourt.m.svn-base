%% PCA

% put your "DATA" is a matrix where each row is a cell and each column a
% parameter
DATA = rand(2000,4); % here I define it to be a random matrix 

[PC, SCORE, LATENT, TSQUARE] = princomp(DATA); %does the PCA

figure(10);clf 
hold off;
set(gcf,'Position',[700 0 700 1120]) % sets position of figure
subplot(2,2,1)
plot(SCORE(:,1),SCORE(:,2),'.')% Scatter plot first two PCs
axis tight;
title('PCA1 vs PCA2');
subplot(2,2,2) % plots amount of variance "explained" by each eigenvector
plot(LATENT,'x'); 
axis tight;
title('EigenValues');
subplot(2,2,3:4)
imagesc(hist3([SCORE(1:end,1) SCORE(1:end,2)],[200 200]));
axis xy; axis tight;
colorbar;colormap(gca,'hot')
title('PCA1 vs PCA2');


%% kmeans
% It is setup to run kmeans in PCA space, but you can use all the principle components so as not to 
% loose any information. Running it in PCA space however allows you to plot
% the first 2  (or 3 PCs if you slightly modify the plotting) for visualization.

nk = 16; % number of clusters
nPC = size(DATA,2); % number of principle components to use (set to number of columns in data to use all data)

[klusters Center]= kmeans(SCORE(1:end,1:nPC),nk); 


% Show klusters in 2dim PCA space
tt = size(colormap,1);temp = colormap;my_colors=temp([1:round(tt/10):tt],:);my_colors = [my_colors; my_colors; my_colors;my_colors;my_colors;my_colors]; % just helpful for displaying different color dots
figure(20);clf;
hold off; set(gcf,'Position',[700 0 1500 1120]);
subplot(2,1,1)
clear sl;
for i =1:nk
    pp=plot(SCORE(klusters==i,1),SCORE(klusters==i,2),'.');
    set(pp,'color',my_colors(i,:));
    hold all
    sl{i} = num2str(i); % used to create legend
    
    %     pause; % UNCOMMENT if there lots of clusters it can be helpful to see them
    %     one at a time with a pause
    
end
legend(sl);
title(['Overclustered in PCA1 vs PCA2 (nPC=' num2str(nPC) ')']);

%Create and dendrogram
Y = pdist(Center,'mahalanobis'); % find the distance between each cluster and all other clusters
                                % compute it in terms of mahalanobis distance
Z = linkage(Y);  % this function is being used to put data in right form for dendrogram has additional options, read help
subplot(2,1,2); % plot dendrogram
dendrogram(Z,0,'orientation','top'); % 
xlabel('cluster #')

% Now you can manually joint clusters based on the PCA projection and the
% dendrogram