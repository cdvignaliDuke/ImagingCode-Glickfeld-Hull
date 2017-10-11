function laser_on_ind_conserv = conservative_laser_on(laser_power_vec_ttl);
%function to generate a much more conservative estimate of when the laser
%on. Only to be used for purposes of extracting signal (PCA/ICA and thresh)
%it will absolutely label some powered frames as unpowered. 

%find the frames where the laser first turns on for each trial
laser_up_ind = zeros(1,length(laser_power_vec_ttl));
for frame_num = 2:length(laser_power_vec_ttl)-1
    if laser_power_vec_ttl(frame_num) == 1
        if laser_power_vec_ttl(frame_num-1)==0
            if laser_power_vec_ttl(frame_num+1)==1
                laser_up_ind(1,frame_num) = 1;
            end
        end
    end
end

laser_up_ind = find(laser_up_ind);
laser_down_ind = nan(1,length(laser_up_ind));
%if laser on
if laser_power_vec_ttl(end)==1
    laser_up_ind = laser_up_ind(1:end-1);
    laser_down_ind = laser_down_ind(1:end-1);
end

for ii = 1:length(laser_up_ind)
    laser_down_ind(ii) = find(laser_power_vec_ttl(laser_up_ind(ii):end)==0 ,1,'first') + laser_up_ind(ii) -1;
end

laser_up_ind = laser_up_ind+20;
laser_down_ind = laser_down_ind-20;

laser_on_ind_conserv = zeros(1,length(laser_power_vec_ttl));
for ii = 1:length(laser_up_ind)
    laser_on_ind_conserv(laser_up_ind(ii):laser_down_ind(ii)) =1;
end
laser_on_ind_conserv = find(laser_on_ind_conserv);

end 





