
filelist='exp_list_scmr8_4.xls';
root='\\zmey\storlab\data\Soumya\scmr8';
overwrite=0;


s=readXLSfilelist(filelist);


mkdir(root, 'stats');
stats_dir = [root,'\stats'];
lowcut_filter=120; % frames  
tcourse_dir=[root, '\tcourse'];


for i=1:length(s)
    
    % if you want to add 'r' in front of filename
    % s(i).filename=['r',s(i).filename,'_ch1'];
    if isempty(findstr(s(i).stim_code,'ori'))
        continue;
    end
    i
    tcoursefname=[s(i).filename,'_tcourse.mat'];
    statsfname=[s(i).filename,'_Oristats.mat'];
    if overwrite==0 && exist(fullfile(stats_dir, statsfname))==2
        continue;
    end
    
    % load tcourses
    load (fullfile(tcourse_dir, tcoursefname));
    
    % setup parameters
    nframes=size(timeCourses,1);
    expt=[];
    expt.dummy=[];
    expt = setupstim(expt,s(i).Non, s(i).Noff, s(i).Nstim_per_run, s(i).Nframes);
    expt.stims = getepochs (s(i).Noff, s(i).Non, s(i).Nstim_per_run, s(i).Nrep,1);

    % preprocessing of time courses
    timeCourses(expt.dur+1:end,:)=[];
    lowcutTimeCourses = tcLowCut (timeCourses, lowcut_filter, 'gaussian', 1);
    DC = repmat(mean(lowcutTimeCourses,1),expt.dur,1);
    normTimeCourses = (lowcutTimeCourses-DC)./DC*100;   % percent signal change
    data_tables = tcEpochAverage(lowcutTimeCourses,expt.stims);

    % get stats
    Oristats = tcOriStats(data_tables, 0.05);
    
    % save stats
    save (fullfile(stats_dir, statsfname) ,'Oristats');
end

