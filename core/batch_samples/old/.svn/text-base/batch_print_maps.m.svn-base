filelist='exp_list_scmr6_2.xls';
s=readXLSfilelist(filelist);

data_dir='\\zmey\storlab\data\Soumya\scmr6';
tif_dir=[data_dir,'\tif'];

map_root_dir = 'E:\users\kenichi\scmr6\pixmaps';
label_dir = 'E:\users\kenichi\scmr6\labelimg'
for i=1:length(s)
    fname=s(i).filename;
    if fname(1)=='r'
        fname(1)=[];
    end
    load(fullfile(tif_dir,[fname,'.mat']));
    [expinfo, stiminfo, scaninfo]=ExpInfo2txt(s(i),params) 
    
    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);

    htxt1=annotation('textbox',[0,0.75,0.33,0.25]);
    htxt2=annotation('textbox',[0.33,0.75,0.33,0.25]);
    htxt3=annotation('textbox',[0.66,0.75,0.33,0.25]);
    % left bottom corner (0,0.5)  box size (0.5,0.45) ;  x-> hotizontal, y-> vertical (upward)

    set(htxt1,'FontSize',12);
    set(htxt2,'FontSize',12);
    set(htxt3,'FontSize',12);
    set(htxt1,'string',expinfo);
    set(htxt2,'string',stiminfo);
    set(htxt3,'string',scaninfo);

    map_dir = fullfile(map_root_dir,s(i).filename);

    mapname=['FOV.tif'];
    map=imread(fullfile(map_dir,mapname));
    subplot(10,8,25);
    imshow(map);
    mapname(strfind(mapname, '_'))='-';
    title(mapname);

    mapname=[s(i).filename,'_shadelabel.tif'];
    [map, cmap]=imread(fullfile(label_dir,mapname));
    sz=size(map);
    map=map(1:sz(1),1:sz(2)/2);
    subplot(10,8,26);
    imshow(map,cmap);

    
    
    for stim=1:s(i).Nstim_per_run
        mapname=['dF_', num2str(stim), '.tif'];
        map=imread(fullfile(map_dir,mapname));
        subplot(10,8,stim+32);
        imshow(map*1.5);
        mapname(strfind(mapname, '_'))='-';
        title(mapname);
    end
    for stim=1:s(i).Nstim_per_run
        mapname=['ratio_', num2str(stim), '.tif'];
        map=imread(fullfile(map_dir,mapname));
        subplot(10,8,stim+56);
        imshow(map*1.5);
        mapname(strfind(mapname, '_'))='-';
        title(mapname);
    end
     print ('-dpsc2', '-append','scmr6_maps_hc.ps');
fclose all
    
end
fclose all
