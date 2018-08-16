for ii = 1: length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    raw_F = load([image_dest, '_raw_F.mat']);
    data_tc = raw_F.data_tc;
    %data_tc is the fluorescence data for this session, and it will be a
    %n*num of frames double, n=number of ROIs.
    behav_output = load([behav_dest '\' days '_behavAnalysis.mat']);
    stay = behav_output.frames_stay_cell;
    for n = 1:size(data_tc,1);
        data_tc_stay = data_tc(n,cell2mat(stay));
        F_staybase = mean(data_tc_stay);
        dfOvF_staybase(n,:) = (data_tc(n,:) - F_staybase)/F_staybase;
    end
   save([image_dest,'_dfOvF_staybase'],'dfOvF_staybase');
end