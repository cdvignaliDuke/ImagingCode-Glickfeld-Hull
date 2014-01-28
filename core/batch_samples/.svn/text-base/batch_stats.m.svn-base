
filelist='exp_list_scmr6_2.xls';
root='\\zmey\storlab\data\Soumya\scmr6';
analysis_dir='\\zmey\storlab\data\Soumya\scmr6';
overwrite=0;

s=readXLSfilelist(filelist);

dirs=GetDirs(data_dir, analysis_dir);

for i=1:length(s)
    % load tcourses
    fnames=GetFnames(s(i).fname,data_dir,analysis_dir);
    load (fnames.tcourse);
    
    % setup parameters
    nframes=size(timeCourses,1);
    expt=[];
    expt.dummy=[];
    expt = setupstim(expt,s(i).Non, s(i).Noff, s(i).nstim_per_run, nframes);
    expt.stims = getepochs (s(i).Noff, s(i).Non, s(i).nstim_per_run, s(i).Nrep,1);

    [ProcTimeCourses,data_tables]=tcProcess(timeCourses,s(i));
    % get stats
    stats = tcStats(data_tables, 0.05);
    
    % save stats
    save (fnames.stats ,'stats','ProcTimeCourses');
end

