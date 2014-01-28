
% filelist='exp_list_scmr6_2.xls';
root='J:\data\mouse090316';
overwrite=0;


%s=readXLSfilelist(filelist);


mkdir(root, 'stats3');
stats_dir = [root,'\stats3'];
tcourse_dir=[root, '\tcourse3'];

s(1).filename='rvstim4';
s(2).filename='rvstim5';
s(3).filename='rvstim6';
s(4).filename='rvstim9';
s(5).filename='rvstim10';
s(6).filename='rvstim11';
s(7).filename='rvstim12';
s(8).filename='rvstim13';
s(9).filename='rvstim14';
s(10).filename='rvstim15';
s(11).filename='rvstim16';
s(12).filename='rvstim17';
s(13).filename='rvstim18';
s(14).filename='rvstim19';
s(15).filename='rvstim20';
s(16).filename='rvstim21';

for i=1:length(s)
    
    % if you want to add 'r' in front of filename
    % s(i).filename=['r',s(i).filename,'_ch1'];
%     if isempty(findstr(s(i).stim_code,'ori'))
%         continue;
%     end
    i
    
    s(i).Non=8;
    s(i).Noff=8;
    s(i).Nstim_per_run=12;
    s(i).Nrep=10;
    
    if i==2
    s(i).Nrep=3;
    end    
    if i==10
    s(i).Nrep=5;
    end    
    if i==11
    s(i).Nrep=5;
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

