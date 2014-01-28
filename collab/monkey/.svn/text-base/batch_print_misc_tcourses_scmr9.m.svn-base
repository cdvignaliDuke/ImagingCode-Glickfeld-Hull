filelist='exp_list_scmr9_2_misc.xls';
s=readXLSfilelist(filelist);


map_root_dir = '\\zmey\storlab\data\Soumya\scmr9\pixelMiscmaps';
data_dir='\\zmey\storlab\data\Soumya\scmr9\';

outfname='E:\users\kenichi\misc_tcourses_scmr9.ps';

tcourse_dir=[data_dir,'tcourse\'];
label_dir=[data_dir,'labelimg\'];
stats_dir=[data_dir,'stats\'];
tif_dir=[data_dir,'tif\'];

p_th=0.01;
R_th=0.05;

Ncell_per_page=20;


for i=1:length(s)
%for i=1
    

    fname=s(i).filename;
    Non=s(i).Non;
    Noff=s(i).Noff;
    Nstim=s(i).Nstim_per_run;
    nframes_per_run=(Non+Noff)*Nstim;
    fnames=filenameManager(fname);

    %%% print experiment information
    org_fname=fname;
    if org_fname(1)=='r'
        org_fname(1)=[];
    end
    load(fullfile(tif_dir,[org_fname,'.mat']));
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

    % show colormap 
    map=imread(fullfile([map_root_dir,'\',s(i).filename],'dir_dF_HLS_hc.tif'));
    subplot('Position',[0 0 1 0.7]); 
    imshow(map);
    print ('-dpsc2', '-append',outfname);

    % load data

    load(fullfile(stats_dir,fnames.stats));
    load(fullfile(tcourse_dir,fnames.tcourse));
    load(fullfile(label_dir,fnames.labelimg));
    
    filter=find([stats.p_value_resp]<p_th & [stats.R_best_dir]>R_th);
    Ncells=length(filter)

    % process timecoueses
    [Nt,Nallcells]=size(timeCourses);
    ProcTimeCourses=tcProcess(timeCourses,s(i));
    time_ticks=[];
    for j=1:Nstim
        time_ticks=[time_ticks,Noff+1+(Noff+Non)*(j-1)];
    end
    
    % print in multiple pages
    
    Npages=ceil(Ncells./Ncell_per_page);
    for p=1:Npages
        hf=figure;
        set(hf,'PaperUnits','normalized');
        set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);
        offset=Ncell_per_page*(p-1);
        
        % print time courses
        for cell=1:Ncell_per_page
            subplot(8,5,cell)
            if cell+offset>Ncells
                break;
            end
            plot(ProcTimeCourses.norm_avg(:,filter(cell+offset)));
            title(['cell ',num2str(filter(cell+offset))]);
            axis([1,nframes_per_run,-10,30])
            set(gca, 'XTick', time_ticks);
            set(gca, 'XGrid','on');
        end
        
        % print ROIs
        
        cell_mask=map/2+128;

        for cell=1:Ncell_per_page
            if cell+offset>Ncells
                break;
            end
            contour=edge(uint8(labelimg==filter(cell+offset)));
            cell_mask(find(contour==1))=255;
            cell_mask(find(contour==1)+s(i).MATsize*s(i).MATsize)=0;
            cell_mask(find(contour==1)+s(i).MATsize*s(i).MATsize*2)=0;

        end

        subplot('Position',[0 0 1 0.5]); 
        imshow(cell_mask);

        xy=getCellCoordinate(labelimg);
        for j=1:Ncell_per_page
            if j+offset>Ncells
                break;
            end
            num=filter(j+offset);
            text(xy(1,num)+5,xy(2,num),num2str(num));
        end
        print ('-dpsc2', '-append',outfname);
    end
    close all  
end
close all
