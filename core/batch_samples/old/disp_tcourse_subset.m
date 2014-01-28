
%orimap
%filter


filelist='exp_list_scmr6_3.xls';
n=25;   % which run you want to display



s=readXLSfilelist(filelist);

im=orimap;

data_dir='\\zmey\storlab\data\Soumya\scmr6\';

tcourse_dir=[data_dir,'tcourse\'];
label_dir=[data_dir,'labelimg\'];
stats_dir=[data_dir,'stats\'];


fname=s(n).filename;
Non=s(n).Non;
Noff=s(n).Noff;
Nstim=s(n).Nstim_per_run;
nframes_per_run=(Non+Noff)*Nstim;

fnames=filenameManager(fname);

load(fullfile(stats_dir,fnames.stats));
%load(fullfile(stats_dir,fnames.Oristats));
load(fullfile(tcourse_dir,fnames.tcourse));
load(fullfile(label_dir,fnames.labelimg));

Ncells=length(filter)

%cell_ind=1:Ncells

cell_ind=1:15
offset=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



[Nt,Nallcells]=size(timeCourses);

ProcTimeCourses=tcProcess(timeCourses,s(n));

time_ticks=[];
for i=1:Nstim
    time_ticks=[time_ticks,Noff+1+(Noff+Non)*(i-1)];
%    time_ticks=[time_ticks,Noff+1+(Noff+Non)*(i-1), (Noff+Non)*i+1];
end

figure;

for i=cell_ind
    subplot(5,3,i)
    if i+offset>Ncells
        break;
    end
    plot(ProcTimeCourses.norm_avg(:,filter(i+offset)));
    title(['cell ',num2str(filter(i+offset))]);
    axis([1,nframes_per_run,-10,30])
    set(gca, 'XTick', time_ticks);
    set(gca, 'XGrid','on');
end

%%%%%%%%%%%%%%%%%%%%%%
% display cell positions

figure;

%subplot(1,2,1);
%imshow(im);

cell_mask=im;

for i=cell_ind
    if i+offset>Ncells
        break;
    end
    contour=edge(uint8(labelimg==filter(i+offset)));
    cell_mask(find(contour==1))=255;
    cell_mask(find(contour==1)+512*512)=0;
    cell_mask(find(contour==1)+512*512*2)=0;

end

%subplot(1,2,2);
imshow(cell_mask);
%colormap(gray)

xy=getCellCoordinate(labelimg);
for i=cell_ind
    if i+offset>Ncells
        break;
    end
    num=filter(i+offset);
    text(xy(1,num)+5,xy(2,num),num2str(num),'Color',[1 1 1]);
end

