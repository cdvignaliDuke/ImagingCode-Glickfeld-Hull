load('\\zmey\storlab\data\Soumya\scmr6\stats\r080112_057_ch1_stats.mat');
load('\\zmey\storlab\data\Soumya\scmr6\labelimg\r080112_057_ch1_labelimg.mat');
%load('\\zmey\storlab\data\Soumya\scmr6\stats\r080112_049_ch1_stats.mat');
%load('\\zmey\storlab\data\Soumya\scmr6\labelimg\r080112_049_ch1_labelimg.mat');

Nstim=10;
Noff=6;
Non=6;

p_th=0.01;
R_th=0.05;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filter=find([stats.R_best_dir]>R_th & [stats.p_value_resp]<p_th);

Ncells=length(filter);

Ncl=min(ceil(Ncells/2),20);    %number of clusters


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing for cluster analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get ratio_change (dF/F) for each stim for each cell

ratio_change=[stats(filter).dir_ratio_change];
ratio_change=reshape(ratio_change,Nstim,Ncells)';

% Get tcourse for each cell

tcourse=ProcTimeCourses.norm_avg;
tcourse=tcourse(:,filter)';
tcourse=(tcourse*2+circshift(tcourse,[0 1])+circshift(tcourse,[0 -1]))/4;

% normalize with peak dF/F or timecourse change

Nt=size(ProcTimeCourses(1).norm_avg,1);
ratio_change_norm=ratio_change./repmat(max(ratio_change,[],2),1,Nstim);
tcourse_norm=tcourse./repmat(max(tcourse,[],2),1,Nt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k-means clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cl=kmeans(ratio_change,Ncl);
cl2=kmeans(tcourse,Ncl);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get average dF/F and tcourse across all cells in each cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avg_ratio_change=zeros(Nstim,Ncl);
avg_tcourse=zeros(Nt,Ncl);

for i=1:Ncl
    cell_ind=find(cl2==i);      % Here clusering by tcourse is used
    if ~isempty(cell_ind)
        avg_ratio_change(:,i)=mean(ratio_change(cell_ind,:),1);
        avg_tcourse(:,i)=mean(tcourse(cell_ind,:),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reorder cluster numbers by hierarchical clustering
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y=pdist(avg_ratio_change');
Z=linkage(Y);
[H,T,perm]=dendrogram(Z);
close

cl2=perm(cl2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot tcourses of cells in each cluster
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
for i=1:Ncl
    subplot(ceil(Ncl/2),2,i);
    cell_ind=find(cl2==i);
    if ~isempty(cell_ind)
        plot(ProcTimeCourses.norm_avg(:,filter(cell_ind)));
%        plot(tcourse(cell_ind,:)');
        set(gca,'XTick',[Noff:(Noff+Non):(Noff+Non)*Nstim-Non],'XGrid','on'); 
        ylim([-5 30]);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get average dF/F and tcourse again across all cells in each cluster
% Get also cell map
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cellmap=zeros(size(labelimg,1),size(labelimg,2));

for i=1:Ncl
    cell_ind=find(cl2==i);
    if ~isempty(cell_ind)
        avg_ratio_change(:,i)=mean(ratio_change(cell_ind,:),1);
        avg_tcourse(:,i)=mean(tcourse(cell_ind,:),1);
        for j=1:length(cell_ind)
            cellmap(find(labelimg==filter(cell_ind(j))))=i;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correlation matrix across clusters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r=corrcoef(avg_tcourse);
figure;
imagesc(r);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cell maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
imagesc(cellmap);
colorbar;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avg ratio change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
for i=1:Ncl
    subplot(ceil(Ncl/2),2,i);
    bar(avg_ratio_change(:,i));
    ylim([-0.1, 0.2]);
end
    
    




