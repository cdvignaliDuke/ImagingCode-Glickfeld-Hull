for ii = 1:length(sessions)
    image_dest = [image_dest_base sessions{ii} '\' sessions{ii}];
    meta_data_dir = [image_dest_base sessions{ii} '\' sessions{ii} '_meta_data'];
    cluster_dir = [image_dest_base sessions{ii} '\' sessions{ii} '_cluster'];
    [avg_img, roi_sz,data_tc] = calculate_data_tc_jin(meta_data_dir, cluster_dir);
    rawF_dir = [image_dest, '_raw_F.mat'];
    save(rawF_dir, 'data_tc', 'roi_sz', 'avg_img');%automatically saves data_tc to bxOutputs
end