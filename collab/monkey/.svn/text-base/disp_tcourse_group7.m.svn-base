filelist='exp_list_scmr6_4.xls';
s=readXLSfilelist('exp_list_scmr6_4.xls');


oriented=[232,330,447,489,501,512,513,535,545,568,570,576,583,611,618,647,670,1355,1490];
unoriented=[16,72,82,135,148,154,178,213,242,254,621,864,1009];



data_dir='\\zmey\storlab\data\Soumya\scmr6\';

tcourse_dir=[data_dir,'tcourse\'];
label_dir=[data_dir,'labelimg\'];

figure


for n=20:29

    fname=s(n).filename;
    Non=s(n).Non;
    Noff=s(n).Noff;
    Nstim=s(n).Nstim_per_run;
    nframes_per_run=(Non+Noff)*Nstim;

    fnames=filenameManager(fname);

    load(fullfile(tcourse_dir,fnames.tcourse));
    load(fullfile(label_dir,fnames.labelimg));


    filter=unoriented(1:5);

    Ncells=length(filter)

    [Nt,Nallcells]=size(timeCourses);

    ProcTimeCourses=tcProcess(timeCourses,s(n));

    time_ticks=[];
    for i=1:Nstim
        time_ticks=[time_ticks,Noff+1+(Noff+Non)*(i-1)];
    %    time_ticks=[time_ticks,Noff+1+(Noff+Non)*(i-1), (Noff+Non)*i+1];
    end


    for i=1:Ncells
        subplot(10,Ncells,i+Ncells*(n-20))
        plot(ProcTimeCourses.norm_avg(:,filter(i)));
        title(['cell ',num2str(filter(i))]);
        axis([1,nframes_per_run,-10,30])
        set(gca, 'XTick', time_ticks);
        set(gca, 'XGrid','on');
    end



end