
filelist='exp_list_scmr8_2.xls';
s=readXLSfilelist(filelist);

data_dir='\\zmey\storlab\data\Soumya\scmr8\';

outfname='cluster_scmr8.ps';
outmatname='cluster_scmr8.mat';

p_th=0.01;
R_th=0.05;

label_dir=[data_dir,'labelimg\'];
stats_dir=[data_dir,'stats\'];

Clusters=cell(length(s),1);

for i=1:length(s)
    Skip=[];
    
    if length(find(Skip==i)) 
        continue
    end
    
    if (isempty(strfind(s(i).stim_code,'Img_flash')) & isempty(strfind(s(i).stim_code,'sweepbar')))
        continue;
    end

    fname=s(i).filename;
    Non=s(i).Non;
    Noff=s(i).Noff;
    Nstim=s(i).Nstim_per_run;
    nframes_per_run=(Non+Noff)*Nstim;

    load(fullfile(stats_dir, [fname,'_stats.mat']));
    load(fullfile(label_dir, [fname,'_labelimg.mat']));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    filter=find([stats.R_best_dir]>R_th & [stats.p_value_resp]<p_th);

    Ncells=length(filter);

    if Ncells<=2

        continue;
    end

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

    for j=1:Ncl
        cell_ind=find(cl2==j);      % Here clusering by tcourse is used
        if ~isempty(cell_ind)
            avg_ratio_change(:,j)=mean(ratio_change(cell_ind,:),1);
            avg_tcourse(:,j)=mean(tcourse(cell_ind,:),1);
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

    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);

    fname2=fname;
    fname2(strfind(fname2,'_'))='-';
    htxt1=annotation('textbox',[0.4,0.95,0.2,0.05]);
    set(htxt1,'FontSize',12);
    set(htxt1,'string',fname2);

    for j=1:Ncl
        subplot(ceil(Ncl/2),2,j);
        cell_ind=find(cl2==j);
        if ~isempty(cell_ind)
            plot(ProcTimeCourses.norm_avg(:,filter(cell_ind)));
    %        plot(tcourse(cell_ind,:)');
            set(gca,'XTick',[Noff:(Noff+Non):(Noff+Non)*Nstim-Non],'XGrid','on'); 
            ylim([-5 30]);
            title(['cluster ',num2str(j)]);
        end
    end

    print ('-dpsc2', '-append',outfname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Get average dF/F and tcourse again across all cells in each cluster
    % Get also cell map
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cellmap=zeros(size(labelimg,1),size(labelimg,2));

    for j=1:Ncl
        cell_ind=find(cl2==j);
        if ~isempty(cell_ind)
            avg_ratio_change(:,j)=mean(ratio_change(cell_ind,:),1);
            avg_tcourse(:,j)=mean(tcourse(cell_ind,:),1);
            for k=1:length(cell_ind)
                cellmap(find(labelimg==filter(cell_ind(k))))=j;
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % disp avg ratio change
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);

    for j=1:Ncl
        subplot(ceil(Ncl/2),2,j);
        bar(avg_ratio_change(:,j));
        ylim([-0.1, 0.2]);
        title(['cluster ',num2str(j)]);
    end

    print ('-dpsc2', '-append',outfname);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Correlation matrix across clusters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);

    r=corrcoef(avg_ratio_change);
    subplot(2,1,1);
    imagesc(r);
    axis square
    colorbar;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cell maps
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(2,1,2);
    imagesc(cellmap);
    axis square
    colorbar;

    print ('-dpsc2', '-append',outfname);
    Clusters{i}.numbers = cl2;
    Clusters{i}.avg_ratio_change=avg_ratio_change;
    Clusters{i}.avg_tcourse=avg_tcourse;
    Clusters{i}.fname=fname;
end    

save(outmatname,'Clusters');



