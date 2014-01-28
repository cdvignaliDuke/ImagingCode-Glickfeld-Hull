
% filelist='exp_list_scmr6_2.xls';
root='E:\users\kenichi\leicatemp\';
overwrite=0;


%s=readXLSfilelist(filelist);


mkdir(root, 'stats');
stats_dir = [root,'\stats'];
tcourse_dir=[root, '\tcourse'];
s(1).filename='rvstim13';
s(2).filename='rvstim14';
s(3).filename='rvstim15';
s(4).filename='rvstim16';
s(5).filename='rvstim17';
s(6).filename='rvstim18';
s(7).filename='rvstim19';
s(8).filename='rvstim20';

for i=1:length(s)
    
    % if you want to add 'r' in front of filename
    % s(i).filename=['r',s(i).filename,'_ch1'];
%     if isempty(findstr(s(i).stim_code,'ori'))
%         continue;
%     end
    i
    
    s(i).Non=10;
    s(i).Noff=10;
    s(i).Nstim_per_run=12;
    s(i).Nrep=10;
    
    if i==6
    s(i).Nrep=7;
    end    
    tcoursefname=[s(i).filename,'_tcourse.mat'];
    statsfname=[s(i).filename,'_Oristats.mat'];
    if overwrite==0 && exist(fullfile(stats_dir, statsfname))==2
        continue;
    end
    
    if exist(fullfile(tcourse_dir, tcoursefname))~=2
        continue;
    end
    
    % load tcourses

    load (fullfile(tcourse_dir, tcoursefname));
    
    [ProcTimeCourses,data_tables]=tcProcess(timeCourses,s(i),30,1);
    % get stats
    Oristats = tcOriStats(data_tables, 0.05);
    
    % save stats
    save (fullfile(stats_dir, statsfname) ,'Oristats');
end

