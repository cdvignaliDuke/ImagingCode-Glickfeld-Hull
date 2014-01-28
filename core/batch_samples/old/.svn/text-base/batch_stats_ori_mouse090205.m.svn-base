
% filelist='exp_list_scmr6_2.xls';
root='\\zmey\storlab\data\Kenichi\mouse090205\';
overwrite=0;


%s=readXLSfilelist(filelist);


mkdir(root, 'stats2');
stats_dir = [root,'\stats2'];
tcourse_dir=[root, '\tcourse'];

s(1).filename='rvstim5';
s(2).filename='rvstim8';
s(3).filename='rvstim12';
s(4).filename='rvstim13';
s(5).filename='rvstim14';
s(6).filename='rvstim15';
s(7).filename='rvstim16';
s(8).filename='rvstim17';
s(9).filename='rvstim18';
s(10).filename='rvstim19';
s(11).filename='rvstim20';
s(12).filename='rvstim21';
s(13).filename='rvstim22';
s(14).filename='rvstim23';
s(15).filename='rvstim24';
s(16).filename='rvstim25';
s(17).filename='rvstim26';

for i=1:length(s)
    
    % if you want to add 'r' in front of filename
    % s(i).filename=['r',s(i).filename,'_ch1'];
%     if isempty(findstr(s(i).stim_code,'ori'))
%         continue;
%     end
    i
    
    s(i).Non=8;
    s(i).Noff=8;
    s(i).Nstim_per_run=8;
    s(i).Nrep=10;
    
    if i==5
    s(i).Nrep=8;
    end    
    
    tcoursefname=[s(i).filename,'_tcourse.mat'];
    statsfname=[s(i).filename,'_Oristats.mat'];
    if overwrite==0 && exist(fullfile(stats_dir, statsfname))==2
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

