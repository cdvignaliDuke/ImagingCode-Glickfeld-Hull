filelist='exp_list_scmr6_2.xls';
s=readXLSfilelist(filelist);

data_dir='\\zmey\storlab\data\Soumya\scmr6';
tif_dir=[data_dir,'\tif'];

map_root_dir = 'E:\users\kenichi\scmr6\pixmaps';
label_dir = 'E:\users\kenichi\scmr6\labelimg'
%for i=1:length(s)
for i=1:1
    fname=s(i).filename;
    if fname(1)=='r'
        fname(1)=[];
    end
    load(fullfile(tif_dir,[fname,'.mat']));
    fileinfo=struct2txt(s(i));
    MPinfo=MPparam2txt(params);
    
    
    hf=figure;
    set(hf,'PaperUnits','normalized');
    set(hf,'PaperPosition',[0.05 0 0.9 1]);

    htxt=annotation('textbox',[0,0.75,1,0.25]);
    % left bottom corner (0,0.5)  box size (0.5,0.45) ;  x-> hotizontal, y-> vertical (upward)

    set(htxt,'string',[fileinfo,MPinfo]);

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
        imshow(map);
        mapname(strfind(mapname, '_'))='-';
        title(mapname);
    end
    for stim=1:s(i).Nstim_per_run
        mapname=['ratio_', num2str(stim), '.tif'];
        map=imread(fullfile(map_dir,mapname));
        subplot(10,8,stim+56);
        imshow(map);
        mapname(strfind(mapname, '_'))='-';
        title(mapname);
    end
     print ('-dpsc2', [s(i).filename,'_maps.ps']);


end
