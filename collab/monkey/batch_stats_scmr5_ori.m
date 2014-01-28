
filelist='exp_list_scmr5_2.xls';
root='\\zmey\storlab\data\Soumya\scmr5';
overwrite=1;


s=readXLSfilelist(filelist);


stats_dir = '\\zmey\storlab\users\monkey_archive\stats\scmr5';
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
    
    [ProcTimeCourses,data_tables]=tcProcess(timeCourses,s(i));
    % get stats
    Oristats = tcOriStats(data_tables, 0.05);
    
    % save stats
    save (fullfile(stats_dir, statsfname) ,'Oristats','ProcTimeCourses');
end

