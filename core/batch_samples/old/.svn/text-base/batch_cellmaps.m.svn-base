
filelist='scmr7colorfilelist.xls';
s=readXLSfilelist(filelist);

mkdir('E:\users\kenichi\scmr7', 'cellmaps');
stats_dir = 'E:\users\kenichi\scmr7\stats';
lowcut_filter=120; % frames  
labelimg_dir='E:\users\kenichi\scmr7\labelimg';
cellmaps_dir='E:\users\kenichi\scmr7\cellmaps';

colormap=[[0.5,0.5,0.5];
    [1,1,1];    %'LonMon'
    [0,1,0];    %'Mon'
    [0,1,0.5];  %'MonLoff'
    [0,1,1];    %'Loff'
    [0,0,0];    %'LoffMoff'
    [1,0,1];    %'Moff'
    [1,0,0.5];  %'LonMoff'
    [1,0,0];    %'Lon'
    [0,0,1];    %'Son'
    [1,1,0]];   %'Soff'

contrast = [20.6, 8.4, 8.5, 8.8, 35.7, 9.7, 8.3, 7.3, 40, 80]/100;

for i=1:length(s)
    % load 
    statsfname=[s(i).fname,'_stats.mat'];
    load (fullfile(stats_dir, statsfname));
    labelimgfname=[s(i).fname,'_labelimg.mat'];
    load (fullfile(labelimg_dir, labelimgfname));
    
    % bestcolor
    
   cellmapfname=[s(i).fname,'_bestcolor_cellmap.tif'];
    cell_to_exclude=find([stats.p_value_resp]>=0.05);
    map=ezCellMap([stats.best_dir],labelimg,cell_to_exclude);
    imwrite(uint8(map),colormap,fullfile(cellmaps_dir, cellmapfname),'tif');

   cellmapfname=[s(i).fname,'_pchange_cellmap.tif'];
    map=ezCellMap([stats.R_best_dir],labelimg);
    jet2=jet;
    jet2(1,:)=[0,0,0];
    map(find(map>0.2))=0.2;
    imwrite(map*5*63,jet2,fullfile(cellmaps_dir, cellmapfname),'tif');

   cellmapfname=[s(i).fname,'_normalized_bestcolor_cellmap.tif'];
   cellmapfname2=[s(i).fname,'_normalized_bestresp_cellmap.tif'];
    cell_to_exclude=find([stats.p_value_resp]>=0.05);
    bestcolor=zeros(length(stats),1);
    bestresp=zeros(length(stats),1);
    for i=1:length(stats)
        norm_resp=stats(i).dir_ratio_change./contrast;
        [resp,best]=max(norm_resp);
        bestcolor(i)=best;
        bestresp(i)=resp;
    end
    map=ezCellMap(bestcolor,labelimg,cell_to_exclude);
    imwrite(uint8(map),colormap,fullfile(cellmaps_dir, cellmapfname),'tif');
    map=ezCellMap(bestresp,labelimg);
    jet2=jet;
    jet2(1,:)=[0,0,0];
    map(find(map>0.2))=0.2;
    imwrite(map*5*63,jet2,fullfile(cellmaps_dir, cellmapfname2),'tif');

end

