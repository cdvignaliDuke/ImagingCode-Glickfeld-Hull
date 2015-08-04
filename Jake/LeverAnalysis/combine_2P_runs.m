clear all
Jake_2P_exptlist
for id = 1:size(date_mat,1)
    date = date_mat(id,:);
    mouse = mouse_mat(id,:);
    subNum = subNum_mat(id,:);
    time = time_mat(id,:);
    nrun = nrun_mat(id,:);
    run = run_mat(id,:,1);
    
    if nrun>1
        disp(['Combining runs: ' date ' ' mouse])
        run_name = [date '_' mouse '_run' run(length(run)-2:end)];
        out_path = fullfile(out_base,run_name);
        dest =  fullfile(out_path,run_name);
        dest_sub = fullfile([dest '_nosub'], [run_name '_nosub']);
        load([dest_sub '_release_movies.mat'])
        load([dest_sub '_press_movies.mat'])
        fail_movie_temp = fail_movie;
        success_movie_temp = success_movie;
        press_movie_temp = press_movie;
        press_long_movie_temp = press_long_movie;
        load([dest_sub '_spont_events.mat']);
        load([dest_sub '_evoked_events.mat'])
        events_temp = events;
        success_temp = success;
        fail_temp = fail;
        press_temp = press;
        data_tc_spont_temp = data_tc_spont;
        fr_lever_temp = fr_lever;
        events_rate_temp = events_rate;
        success_tc_temp = success_tc;
        fail_tc_temp = fail_tc;
        press_tc_temp = press_tc;
        for i = 2:nrun
            run = run_mat(id,:,i);
            run_name = [date '_' mouse '_run' run(length(run)-2:end)];
            out_path = fullfile(out_base,run_name);
            dest =  fullfile(out_path,run_name);
            dest_sub = fullfile([dest '_nosub'], [run_name '_nosub']);
            load([dest_sub '_release_movies.mat'])
            load([dest_sub '_press_movies.mat'])
            fail_movie_temp = cat(1,fail_movie_temp,fail_movie);
            success_movie_temp = cat(1,success_movie_temp,success_movie);
            press_movie_temp = cat(1,press_movie_temp,press_movie);
            press_long_movie_temp = cat(1,press_long_movie_temp,press_long_movie);
            load([dest_sub '_spont_events.mat']);
            load([dest_sub '_evoked_events.mat'])
            events_temp = [events_temp; events];
            success_temp = [success_temp; success];
            fail_temp = [fail_temp; fail];
            press_temp = [press_temp; press];
            data_tc_spont_temp = [data_tc_spont_temp; data_tc_spont];
            fr_lever_temp = [fr_lever_temp fr_lever];
            events_rate_temp = [events_rate_temp; events_rate];
            success_tc_temp = cat(3, success_tc_temp, success_tc);
            fail_tc_temp = cat(3, fail_tc_temp, fail_tc);
            press_tc_temp = cat(3, press_tc_temp, press_tc);
        end
        press_movie = press_movie_temp;
        press_long_movie =  press_long_movie_temp;
        success_movie = success_movie_temp;
        fail_movie = fail_movie_temp;
        fr_lever = fr_lever_temp;
        data_tc_spont = data_tc_spont_temp;
        success_tc = success_tc_temp;
        fail_tc = fail_tc_temp;
        press_tc = press_tc_temp;
        events_rate = mean(events_rate_temp,1);
        events = concatenateStructures(events_temp);
        success = concatenateStructures(success_temp);
        fail = concatenateStructures(fail_temp);
        press = concatenateStructures(press_temp);
        run = run_mat(id,:,1);
        run_name = [date '_' mouse '_run' run(length(run)-2:end) '-00' num2str(nrun-1)];
        out_path = fullfile(out_base,run_name);
        dest = fullfile(out_path,run_name);
        dest_sub = fullfile([dest '_nosub'], [run_name '_nosub']);
        save([dest_sub '_press_movies.mat'], 'press_long_movie', 'press_movie', 'ifi', 'pre_release_frames','post_release_frames');
        save([dest_sub '_release_movies.mat'], 'success_movie', 'fail_movie', 'pre_release_frames','post_release_frames','ifi');
        save([dest_sub '_spont_events'], 'events', 'data_tc_spont', 'fr_lever', 'thresh', 'events_rate');
        save([dest_sub '_evoked_events.mat'], 'success', 'fail', 'press', 'success_tc', 'fail_tc', 'press_tc')
    end
end
