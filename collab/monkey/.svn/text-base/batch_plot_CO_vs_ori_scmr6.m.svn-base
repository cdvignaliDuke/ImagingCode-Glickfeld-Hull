filelist='exp_list_scmr6_3.xls';
overwrite=0;


s=readXLSfilelist(filelist);


stats_dir = '\\zmey\storlab\users\monkey_archive\stats\scmr6';
label_dir = '\\zmey\storlab\data\Soumya\scmr6\labelimg'
CO_dir='\\zmey\storlab\users\monkey_archive\CO analysis files\scmr6\canonical';
COname='scmr6_FOV';

outfname='plot_CO_vs_ori_scmr6.ps';

for i=1:length(s)
    
    % if you want to add 'r' in front of filename
    % s(i).filename=['r',s(i).filename,'_ch1'];
    if isempty(findstr(s(i).stim_code,'ori8')) & isempty(findstr(s(i).stim_code,'ori16')) 
        continue;
    end
    if isnan(s(i).group)
        continue;
    end
    
    % load stats
    statsfname=[s(i).filename,'_Oristats.mat'];
    load (fullfile(stats_dir, statsfname));

    % load labelimg
    labelfname=[s(i).filename,'_labelimg.mat'];
    load (fullfile(label_dir, labelfname));
    
    % load CO
    COfname = [COname, num2str(s(i).group), '.tif']; 
    blobfield=imread (fullfile(CO_dir, COfname));
    
    COvalue=getBlobIndex(blobfield,labelimg);

    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);
    
    subplot(3,1,1);
    COmap=ezCellMap(COvalue, labelimg,[],0);
    imagesc(COmap);
    axis image;
    colorbar;
    title([s(i).filename, ' CO index']);

    subplot(3,1,2);
    ori_width=[Oristats.ori_tuning_width];
    ori_width(find(isnan(ori_width)))=90;
    tw_map=ezCellMap(ori_width, labelimg,[],0);
    imagesc(tw_map);
    axis image;
    colorbar;
    title('ori tuning width (deg)');

    subplot(3,1,3);
    ori_angle=[Oristats.ori_vector_angle];
    ori_map=ezCellMap(ori_angle, labelimg,[],0);
    imagesc(ori_map);
    axis image;
    colorbar;
    title('orientation map (deg)');

    
    print ('-dpsc2', '-append',outfname);
    close all
    
    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);
    
    plot_CO_vs_ori(Oristats,COvalue);
    print ('-dpsc2', '-append',outfname);
    close all

    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);
    
    hist_CO_vs_ori(Oristats,COvalue);
    print ('-dpsc2', '-append',outfname);
    close all

end

