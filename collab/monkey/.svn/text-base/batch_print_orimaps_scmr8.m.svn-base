filelist='exp_list_scmr8_2.xls';
s=readXLSfilelist(filelist);

map_root_dir = '\\zmey\storlab\data\Soumya\scmr8\pixelOrimaps';
data_dir='\\zmey\storlab\data\Soumya\scmr8\';

outfname='orimaps_scmr8.ps';

tcourse_dir=[data_dir,'tcourse\'];
label_dir=[data_dir,'labelimg\'];
stats_dir=[data_dir,'stats\'];
tif_dir=[data_dir,'tif\'];

p_th=0.01;
R_th=0.05;

Nmap_per_page=6;

hf=figure;
set(hf,'PaperUnits','normalized');
set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);

n=0;

for i=1:length(s)

    if isempty(strfind(s(i).stim_code,'Img_ori'))
        continue;
    end
    n=n+1;

    fname=s(i).filename;
    Non=s(i).Non;
    Noff=s(i).Noff;
    Nstim=s(i).Nstim_per_run;
    nframes_per_run=(Non+Noff)*Nstim;
    fnames=filenameManager(fname);

    %%% print experiment information

    txt=[regexprep(s(i).filename,'_','-'), ',', s(i).stim_description, ',depth=' num2str(s(i).depth),',',num2str(s(i).zoom),'x'];

    subplot(3,2,n)
    % show colormap     
    map=imread(fullfile([map_root_dir,'\',s(i).filename],'ori_dF_HLS_hc.tif'));
    imshow(map);
    title(txt,'FontSize', 8);
    if n==6
        print ('-dpsc2', '-append',outfname);
        n=0;
        hf=figure;
        set(hf,'PaperUnits','normalized');
        set(hf,'PaperPosition',[0.05 0.05 0.9 0.9]);
    end
end

if n>0
    print ('-dpsc2', '-append',outfname);
end
close all  

