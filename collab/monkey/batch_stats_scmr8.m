
filelist='exp_list_scmr8_2.xls';
root='\\zmey\storlab\data\Soumya\scmr8';
overwrite=0;

s=readXLSfilelist(filelist);


mkdir(root, 'stats');
stats_dir = [root,'\stats'];
tcourse_dir=[root, '\tcourse'];


for i=1:length(s)
    
    % if you want to add 'r' in front of filename
    % s(i).filename=['r',s(i).filename,'_ch1'];
    
    
    tcoursefname=[s(i).filename,'_tcourse.mat'];
    statsfname=[s(i).filename,'_stats.mat'];
    if overwrite==0 && exist(fullfile(stats_dir, statsfname))==2
        continue;
    end
    
    % load tcourses
    load (fullfile(tcourse_dir, tcoursefname));
    
    [ProcTimeCourses,data_tables]=tcProcess(timeCourses,s(i));
    % get stats
    stats = tcStats(data_tables, 0.05);
    
    % save stats
    save (fullfile(stats_dir, statsfname) ,'stats','ProcTimeCourses');
end

