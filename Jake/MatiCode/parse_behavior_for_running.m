

function [running, airpuff] = parse_behavior_for_running(b_data)
% calculate  vector of the running at each frame
% get the times of the air puff
% Need to debug andrews code and smooth out the interface between AndrewCode and MatiCode 

ITI_COUNTER =10;
SP_COUNTER =10;

running = [];
airpuff  = [];

% how many frames does each counter represent
iti_counter_to_frame =  b_data.nScansOn/ITI_COUNTER;
sp_counter_to_frame =  b_data.nScansOn/SP_COUNTER;
isi_frames =  b_data.preSoundPauseNFrames + b_data.postSoundPauseNFrames;
% Go over all counters and get the running data get spCopunter and ITI  counter
sp_mat  =[];
iti_mat = [];
for i=1:SP_COUNTER
    name = ['itiCounter' num2str(i)];
    cvec = b_data.(name);
    for j=1:length(cvec) 
        if(isempty(cvec{j})) % FIX ME, check why this happened!!!
            cvec{j} = int64(0);
            warning('check this!!!!');
        end
    end
    c_iti= double(cell2mat(cvec));   

    if(i==1)
        c_iti(1) = NaN;  % The first iti could contain  data before the camera is one
    end
    
    % --- each counter - the data in each counter is sacled by number of frames
    c_iti = repmat(c_iti, iti_counter_to_frame,1)/iti_counter_to_frame;
    %     if(i==10)
    %         % this is due to a bug in the xml code
    %         iti10 = double(cell2mat(b_data.tITIWheelCounter));
    %         c_iti = repmat(iti10, iti_counter_to_frame,1)/iti_counter_to_frame;
    %        % c_iti(:,:) = NaN;
    %     end
    iti_mat(end+1:end+iti_counter_to_frame,:) =  c_iti;
    name= ['spCounter' num2str(i)];
    c_sp = double(cell2mat(b_data.(name)));
    c_sp = repmat(c_sp, sp_counter_to_frame,1)/sp_counter_to_frame;
    sp_mat(end+1:end+sp_counter_to_frame,:) = c_sp;
end


isi_inter_count = double(cell2mat(b_data.tISIWheelCounter));
isi_inter_count = repmat(isi_inter_count, isi_frames,1)/isi_frames;
inter_count(:,:) = isi_inter_count;

running_mat= [iti_mat; inter_count; sp_mat]; % # of rows = # of frames per trial       # of columns = # of trials 


airpuff_mat=  zeros(size(running_mat));
airpuff_mat(end-size(sp_mat,1)+1, :) =1; %  the airpuff is delivered in the begining of the frame

tCounter= double(cell2mat(b_data.tCounter));
running =[];
airpuff =[];
base_len = size(running_mat,1);   % = # of frames per trial
for i=1:size(running_mat,2)
    
    pre_fix =[];
    
    if(i==1)
         % maybe we need to add one or two NaN in this case i.e. pre_fix = [NaN NaN]?
    else
        t_len = tCounter(i) - tCounter(i-1);
        % -- add a null frame when the trials are too long. 
      %  pre_fix = nan(t_len - base_len,1);
        
%         if(t_len>base_len)
%             pre_fix = NaN;
%         end
    end
    vec_run = running_mat(:,i);        % if there were too many pulses associated for a given trial (lumped frames) this code simply adds some NaNs to the beginning of the trial. It does NOT seek out wehre within the trial they are needed. 
    vec_airpuff = airpuff_mat(:,i);
    running = [running; pre_fix; vec_run ];
    airpuff = [airpuff; pre_fix;vec_airpuff ];
end
% running = running_mat(:);
%airpuff = airpuff_mat(:);

return;


