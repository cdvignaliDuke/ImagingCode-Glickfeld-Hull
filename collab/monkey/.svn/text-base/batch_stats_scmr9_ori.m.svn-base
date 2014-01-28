
filelist='exp_list_scmr9_2.xls';
root='\\zmey\storlab\data\Soumya\scmr9';
overwrite=1;


s=readXLSfilelist(filelist);


stats_dir = '\\zmey\storlab\users\monkey_archive\stats\scmr9';
tcourse_dir=[root, '\tcourse'];


for i=1:length(s)
    
    % if you want to add 'r' in front of filename
    % s(i).filename=['r',s(i).filename,'_ch1'];
    if isempty(findstr(s(i).stim_code,'ori8')) & isempty(findstr(s(i).stim_code,'ori16'))
        continue;
    end

    
    
    tcoursefname=[s(i).filename,'_tcourse.mat'];
    statsfname=[s(i).filename,'_Oristats.mat'];
    if overwrite==0 && exist(fullfile(stats_dir, statsfname))==2
        continue;
    end
    
    % load tcourses
    load (fullfile(tcourse_dir, tcoursefname));
    
    [nframes, ncells]=size(timeCourses);     
    expt=[];
    expt.dummy=[];
    expt = setupstim(expt,s(i).Non, s(i).Noff, s(i).Nstim_per_run, nframes);

    
    [ProcTimeCourses,data_tables]=tcProcess(timeCourses,s(i),30,expt.trialdur);

    % get stats
    Oristats = tcOriStats(data_tables, 0.05);
    
    % save stats
    save (fullfile(stats_dir, statsfname) ,'Oristats','ProcTimeCourses');
end

