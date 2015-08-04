function [df_trig_movie, df_f_trig_movie]=  adapt_baselines(trig_movie, ...
    pre_event,use_external, ext_avg) 
% calculate DF  and DF/F from a triggered movie
df_trig_movie = nan(size(trig_movie));
df_f_trig_movie = nan(size(trig_movie));


for i=1:size(trig_movie,1)
    if(exist('use_external','var')&& use_external)
        base_avg = ext_avg;
    else        
        use_inx = 1:pre_event;
        base_avg = mean(trig_movie(i,:,use_inx),3);
    end
    for j=1:size(trig_movie(i,:,:),3)
        df_trig_movie(i,:,j) = trig_movie(i,:,j) - base_avg;
        df_f_trig_movie(i,:,j) = df_trig_movie(i,:,j)./base_avg;
        
    end
end









