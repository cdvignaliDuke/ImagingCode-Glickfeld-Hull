% once time use script for loading day 1 and day N for each dataset, then
% motion registering and taking an average image. then compare the two on a
% pdf
clear;
out_dir = 'Z:\Analysis\Cue_reward_pairing_analysis\2P\FoV_matching_outputs_thresholding\motion_corrected_avg_FoVs\';
%day 1 of 500ms: 
days_1 = {'170416_img90', '170417_img91', '170420_img92', '170510_img93', '170513_img89', '170524_img94', ...
    '170911_img036', '170912_img039', '170921_img044', '171113_img050', '171203_img053', '171227_img067', ...
    '180104_img070', '180322_img077', '180507_img081', '180425_img084', '180509_img085'};

days_1_stable = {[1:200], [1:200], [1:200], [1:200], [1:200], [1:300], ...
    [1:300], [1:300], [1:300], [1:300], [63:320], [1:200], ...
    [1:200], [125:270], [1:200], [1:200], [30:280]};

%post learning 500ms
days_post = {'170426_img90',  '170425_img91', '170428_img92', '170518_img93', '170519_img89', '170529_img94', ...
    '170926_img036', '170917_img039', '170926_img044', '171122_img050', '171207_img053', '180104_img067', ...
    '180108_img070', '180403_img077', '180518_img081', '180502_img084', '180517_img085'};

days_post_stable = {[1:200],  [150:350], [50:250], [1:200], [50:250], [1:300], ...
    [1:300], [335:635], [325:525], [1:200], [1:200], [1:200], ...
    [1:200], [1:200], [600:750], [400:600], [1:300]};

% subjNum_cell = {'990',  '991', '992', '993', '989', '994', ...
%     '036', '039', '044', '050', '053', '067', ...
%     '070', '077', '081', '084', '085'};

subjNum_cell = {'90',  '91', '92', '93', '89', '94', ...
    '036', '039', '044', '050', '053', '067', ...
    '070', '077', '081', '084', '085'};

color_code_qc = {'g', 'y', 'y', 'y', 'r', 'g', 'g', 'r', 'y', 'y', 'r', 'y', 'g', 'y', 'g', 'g', 'y'};  %left of after img067

image_file_num = {'_000', '_001', '_002'};

for session_num = 11 % 1:length(days_1)
    figure;
    
    %% load data
    for pre_post = 1:2
        if pre_post ==1
            this_session = days_1{session_num};
        else
            this_session = days_post{session_num};
        end
        session_date = this_session(1:6);
        subjNum = subjNum_cell{session_num};
        
        %identify sessions with low snr
%         if sum( strcmp(this_session,  {'170911_img036', '170911_img036', '170912_img039', '170917_img039', '170926_img044', '171113_img050',  '171122_img050', ...
%                 '171203_img053', '171207_img053',  '171227_img067', '180104_img067', '180104_img070', '180509_img085', })) == 1;
%         end
        
        for file_order = 1:length(image_file_num)
            fname = ['img', subjNum, '_000', image_file_num{file_order}];
            if exist(['Z:\Data\2P_imaging\', session_date, '_img', subjNum '\img', subjNum, '\', fname, '.mat']) ~= 0
                cd(['Z:\Data\2P_imaging\', session_date, '_img', subjNum '\img', subjNum]);
                %fname_dir = dir(['Z:\Data\2P_imaging\', session_date, '_img', subjNum '\img', subjNum])
                break
            end
        end
        [data, ~, ~] = loadFile(['Z:\Data\2P_imaging\', session_date, '_img', subjNum '\img', subjNum, '\'], image_file_num{file_order}, 1, 1500);
        %data = squeeze(sbxread(fname,0,1500));
        
        %% motion correct
        
%         %determine frame nums to use for selecting stable frame
%         sz = size(data);
%         laser_on_ind = [1:sz(3)];
%         frame_nums_for_ref30 = laser_on_ind(randi([1,size(laser_on_ind,2)],1,30));
%         frame_nums_for_samp100 = laser_on_ind(round(linspace(1,size(laser_on_ind,2))));
%         
%         % select 30 random frames from throughout the movie
%         ref30 = data(:,:,frame_nums_for_ref30);
%         
%         %motion register each of the 30 random frames to 100 frames from the movie. Find the one with the lowerst dshift
%         samp100 = data(:,:,frame_nums_for_samp100);
%         dshift = [];
%         for r = 1:size(ref30,3)
%             [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
%             dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
%         end
%         
%         %pick the frame which had the lowest dshift and motion register full movie to that frame
%         min_f = find(dshift == min(dshift));
%         if length(min_f) > 1
%             min_f = min_f(1,1);
%         end
%         if low_snr == 1  
%             sum_ref30 = sum(sum(ref30,2),1);
%             img_ref_ind = find( sum_ref30 == max(sum_ref30));
%             img_ref = ref30(:,:,img_ref_ind);
%         else
%             img_ref = ref30(:,:,min_f);
%         end
        
        if pre_post == 1
            img_ref = mean(data(:,:,days_1_stable{session_num}),3);
            [reg_out, img_reg] = stackRegister(data, img_ref);
        else
            img_ref = mean(data(:,:,days_post_stable{session_num}),3);
            [reg_out_2, img_reg_2] = stackRegister(data, img_ref);
        end
        
        %% save to cdr
        figure;
        if strcmp(color_code_qc{session_num}, 'r');
            suptitle('RED');
        elseif strcmp(color_code_qc{session_num}, 'y');
            suptitle('YELLOW');
        elseif strcmp(color_code_qc{session_num}, 'g');
            suptitle('GREEN');
        end
        if pre_post == 1
            %subplot(2,1,1);
            img_reg_avg = mean(img_reg,3);
            imagesc(img_reg_avg); truesize;
            title(['day 1: ', subjNum_cell{session_num}]);
            savefig([out_dir, subjNum_cell{session_num}]);
            writetiff(img_reg_avg, [out_dir, subjNum_cell{session_num}]);
        else
            %subplot(2,1,2)
            img_reg_avg = mean(img_reg_2,3);
            imagesc(img_reg_avg); truesize;
            title(['post-learning: ', subjNum_cell{session_num}]);
            savefig([out_dir, subjNum_cell{session_num}, '_2']);
            writetiff(img_reg_avg, [out_dir, subjNum_cell{session_num}, '_2']);
        end
       
    end
    
    %print([out_dir, subjNum_cell{session_num}, '.pdf'], '-dpdf');
end





