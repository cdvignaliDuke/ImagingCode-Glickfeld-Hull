% % get_all_manuscript_data
% % collect all data from WF and 2P expts for the 2019 eLife paper
% 
clear;
list_of_failed_copies = [];
bx_dir = ['Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\'];
WF_img_dir =  ['Y:\home\jake\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
dest_folder = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\'];
dest_folder_WF = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\wide-field data'];

%% find and upload all WF imaging days data
% WF_CRP_list_of_days;
% for day_num = 1:3
%    if day_num ==1
%        curr_dataset = days_1;
%        curr_roi = days_1_ROI;
%        dest_id = 'day_1';
%    elseif day_num==2
%        curr_dataset = days_post;
%        curr_roi = days_post_ROI;
%        dest_id = 'day_N1';
%    elseif day_num ==3
%        curr_dataset = days_UR;
%        curr_roi = days_UR_ROI;
%        dest_id = 'day_N2';
%    end
%    
%    for sess_num = 1:length(curr_dataset)
%        this_session = curr_dataset{sess_num};
%        if length(this_session) <6
%            continue
%        end
%         bx_data = get_bx_data(bx_dir, this_session);
%         load([WF_img_dir, this_session, '_bx_outputs'], 'tc_dfoverf');
%         tc_dfoverf = tc_dfoverf(curr_roi{sess_num},:);
%         save([dest_folder_WF, '\', dest_id, '\', this_session, '_timecourses'], 'tc_dfoverf');
%         save([dest_folder_WF, '\', dest_id, '\', this_session, '_behavior'], 'bx_data');
%    end
% end
% 
% days_list{1} = {'170321_img86', '170322_img86', '170323_img86', '170325_img86', '170326_img86', '170327_img86', '170330_img86'};
% days_list{2} = {'170408_img88', '170410_img88', '170411_img88', '170412_img88', '170413_img88', '170414_img88', '170415_img88', '170417_img88', '170418_img88', '170419_img88', '170420_img88', '170421_img88', '170422_img88', '170423_img88', '170424_img88', '170425_img88', '170426_img88'};
% days_list{3} = {'170605_img95', '170606_img95', '170607_img95', '170608_img95', '170609_img95', '170610_img95', '170611_img95', '170612_img95', '170613_img95', '170614_img95'};
% days_list{4} = {'170612_img96', '170613_img96', '170614_img96', '170615_img96', '170620_img96', '170621_img96', '170622_img96', '170623_img96', '170624_img96'};
% days_list{5} = {'170628_img98', '170629_img98', '170701_img98', '170703_img98', '170704_img98', '170705_img98', '170706_img98'};
% days_list{6} = {'170705_img99', '170706_img99', '170707_img99', '170708_img99', '170710_img99', '170711_img99', '170712_img99'};
% days_to_PL_imaging_session = [4, 7, 10, 8, 5, 5];
% 
% %find and collect all WF training days data
% dest_id = 'training_days';
% for animal_num = 1:length(days_list)
%     this_animal = days_list{animal_num}(1);
%     mouse_num = ['i9', this_animal{1}(end-1:end)];
%     if ~exist([dest_folder_WF, '\', dest_id, '\', mouse_num,  '_behavior'])
%         mkdir([dest_folder_WF, '\', dest_id, '\', mouse_num,  '_behavior']);
%     end
%     
%     for this_sess_num = 1:length(days_list{animal_num})
%         this_session = days_list{animal_num}(this_sess_num);
%         this_session = this_session{1};
%         this_bx_data = get_bx_data(bx_dir, this_session);
%         save([dest_folder_WF, '\', dest_id, '\', mouse_num, '_behavior' '\', this_session], 'this_bx_data');
%         if this_sess_num == days_to_PL_imaging_session(animal_num)
%             break
%         end
%     end
% end

%% need a sections which can collect all 2P training days 
% clear days_list;
% dest_id = 'training_days';
% days_list{1} = {'170513_img89', '170515_img89', '170516_img89', '170517_img89', '170518_img89', '170519_img89', '170522_img89', '170523_img89', '170524_img89', '170525_img89', '170526_img89', '170527_img89', '170529_img89', '170530_img89', '170531_img89'};
% days_list{2} = {'170416_img90', '170417_img90', '170418_img90', '170419_img90', '170422_img90', '170423_img90', '170424_img90', '170425_img90', '170426_img90', '170427_img90', '170428_img90', '170429_img90', '170501_img90', '170502_img90', '170503_img90', '170504_img90', '170505_img90', '170506_img90', '170508_img90', '170509_img90', '170510_img90', '170511_img90', '170512_img90', '170513_img90', '170515_img90'};
% days_list{3} = {'170417_img91', '170418_img91', '170419_img91', '170420_img91', '170422_img91', '170423_img91', '170424_img91', '170425_img91', '170426_img91', '170427_img91', '170428_img91', '170429_img91', '170501_img91', '170502_img91', '170503_img91', '170504_img91'};
% days_list{4} = {'170420_img92', '170422_img92', '170423_img92', '170424_img92', '170425_img92', '170426_img92', '170427_img92', '170428_img92', '170429_img92', '170501_img92', '170503_img92', '170504_img92', '170505_img92', '170506_img92', '170509_img92'};
% days_list{5} = {'170510_img93', '170511_img93', '170512_img93', '170513_img93', '170515_img93', '170516_img93', '170517_img93', '170518_img93', '170519_img93', '170520_img93', '170522_img93', '170523_img93', '170524_img93', '170525_img93'};
% days_list{6} = {'170524_img94', '170525_img94', '170526_img94', '170527_img94', '170529_img94', '170530_img94', '170531_img94', '170601_img94', '170602_img94', '170604_img94', '170605_img94', '170606_img94'};
% days_list{7} = {'170921_img044', '170922_img044', '170923_img044', '170924_img044', '170925_img044', '170926_img044'};
% days_list{8} = {'171113_img050', '171115_img050', '171116_img050', '171119_img050', '171120_img050', '171121_img050', '171122_img050'};
% days_list{9} = {'171227_img067', '171228_img067', '171229_img067', '171230_img067', '180101_img067', '180102_img067', '180103_img067', '180104_img067'};
% days_list{10} = {'180104_img070', '180105_img070', '180106_img070', '180107_img070', '180108_img070'};
% days_list{11} = {'180322_img077', '180323_img077', '180324_img077', '180325_img077', '180326_img077', '180327_img077', '180328_img077', '180329_img077', '180330_img077', '180331_img077', '180401_img077', '180402_img077', '180403_img077', '180404_img077'};
% days_list{12} = {'180507_img081', '180508_img081', '180509_img081', '180510_img081', '180511_img081', '180512_img081', '180513_img081', '180514_img081', '180515_img081', '180516_img081', '180517_img081', '180518_img081', '180521_img081'};
% days_list{13} = {'180425_img084', '180426_img084', '180427_img084', '180428_img084', '180429_img084', '180430_img084', '180501_img084', '180502_img084'};
% days_list{14} = {'180509_img085', '180510_img085', '180511_img085', '180512_img085', '180513_img085', '180514_img085', '180515_img085', '180516_img085', '180517_img085', '180518_img085'};
% days_list{15} = {'181213_img087', '181214_img087', '181215_img087', '181216_img087', '181217_img087', '181218_img087'};
% days_list{16} = {'181214_img088', '181215_img088', '181216_img088', '181217_img088', '181218_img088', '181219_img088', '181220_img088', '181221_img088'};
% days_list{17} = {'181214_img089', '181215_img089', '181216_img089', '181217_img089', '181218_img089', '181219_img089'};
% days_list{18} = {'181213_img091', '181214_img091', '181215_img091', '181216_img091', '181217_img091', '181218_img091', '181219_img091', '181220_img091'};
% days_list{19} = {'190213_img092', '190214_img092', '190215_img092', '190216_img092', '190217_img092', '190218_img092', '190219_img092', '190220_img092', '190221_img092', '190222_img092', '190223_img092', '190224_img092', '190225_img092', '190226_img092'};
% days_list{20} = {'190215_img093', '190216_img093', '190217_img093', '190218_img093', '190219_img093', '190220_img093', '190222_img093', '190225_img093', '190226_img093', '190227_img093', '190228_img093', '190301_img093', '190302_img093', '190304_img093', '190306_img093', '190308_img093'};
% days_list{21} = {'190310_img094', '190311_img094', '190312_img094', '190313_img094', '190314_img094', '190315_img094', '190316_img094', '190317_img094', '190318_img094', '190319_img094', '190320_img094', '190321_img094', '190325_img094', '190327_img094'};
% 
% days_to_PL_imaging_session = [6, 11, 8, 8, 8, 5, ...
%     6, 7, 8, 4, 13, 12, 8, ...
%     9, 5, 8, 5, 5, 5, 5, 6];
% all_bx_animals = {'i989', 'i990', 'i991', 'i992', 'i993', 'i994', ...
%     'i044', 'i050', 'i067', 'i070', 'i077', 'i081', 'i084', ...
%     'i085', 'i087', 'i088', 'i089', 'i091', 'i092', 'i093', 'i094'};
% 
% for this_lobule = 1:3
%     %Select the appropriate folder 
%     if this_lobule ==1
%         dest_folder_2P = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\two-photon LS\'];
%         animal_subset = [];
%         CC_paper_expts_LS;
%     elseif this_lobule==2
%         dest_folder_2P = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\two-photon C1\'];
%         animal_subset = [];
%         CC_paper_expts_Crus;
%     elseif this_lobule==3 
%         dest_folder_2P = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\two-photon C2\'];
%         animal_subset = [];
%         CC_paper_expts_Crus;
%     end
%     
%     %find all the training days for each animal
%     curr_dataset = expt(1);
%     for animal_num = days_list_subset; %1:length(days_list)
%         this_animal_start = days_list{animal_num}(1);
%         mouse_num = [all_bx_animals{animal_num}];
%         if ~exist([dest_folder_WF, '\', dest_id, '\', mouse_num, '_behavior'])
%             mkdir([dest_folder_WF,  '\', dest_id, '\', mouse_num, '_behavior']);
%         end
%         
%         for this_sess_num = 1:length(days_list{animal_num})
%             this_session = days_list{animal_num}(this_sess_num);
%             this_bx_data = get_bx_data(bx_dir, this_session);
%             save([dest_folder_WF,  '\', dest_id, '\', mouse_num, '\', this_session, '_behavior'], 'this_bx_data');
%             if this_session_num == days_to_PL_imaging_session(animal_num)
%                 break
%             end
%         end
%     end
% end

%% 2P data LS sessions

% CC_paper_expts_LS;  
% TP_img_dir = ['Y:\home\lindsey\Analysis\2P\Jake\'];
% bx_dir = ['Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\'];
% dest_folder_2P = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\two-photon LS\'];
% load([dest_folder, 'expNums_2P']);
% lobule_ID = 3;
% 
% for day_num = 1:3
%     %determine the behavior condition of the imaging session
%    if day_num ==1
%        curr_dataset = expt(day_num);
%        dest_id = 'day_1';
%    elseif day_num==2
%        curr_dataset = expt(day_num);
%        dest_id = 'day_N1';
%    elseif day_num ==3
%        curr_dataset = expt(day_num);
%        dest_id = 'day_N2';
%    end
%    sess_included =  expNums(day_num).mouse_name([~isnan(expNums(day_num).numDendrites(:,lobule_ID))],1); 
%    
%    %for each session included in that behavior condition...
%    for sess_num = 1:length(curr_dataset.mouse)
%        %select mouse
%        mouse = strtrim(curr_dataset.mouse(sess_num,:));
%        date = curr_dataset.date(sess_num,:);
%        this_session = [date, '_', mouse];
%        
%        %determine whether this mouse is included in this bx condition
%        mouse_rep = cell(size(sess_included,1),1);
%        mouse_rep(:) = {mouse};
%        mouse_ind = find(cellfun(@strcmp, sess_included, mouse_rep));
%        if isempty(mouse_ind)
%            continue
%        end
%        
%        %load mask and bx file
%        img_fn2 = [date '_img' mouse];
%        bx_data = get_bx_data(bx_dir, img_fn2);
%        load([TP_img_dir, this_session, '\', this_session, '_ROI_TCs'], 'mask3D');
%        load([TP_img_dir, this_session, '\', this_session, '_targetAlign.mat']);
%        
%        %copy and save data
%        save([dest_folder_2P,  '\', dest_id, '\', this_session, '_mask3D'], 'mask3D');
%        save([dest_folder_2P, '\',  dest_id, '\', this_session, '_behavior'], 'bx_data');
%        save([dest_folder_2P, '\', dest_id, '\', this_session, '_time_courses'], 'targetAligndFoverF', 'ind_omit', ...
%            'ind_omit_long', 'ind_omit_short', 'ind_rew', 'ind_rew_postomit', 'ind_rew_preomit', 'ind_unexp', 'ind_unexp_long', 'ind_unexp_short', ...
%            'postwin_frames', 'prewin_frames', 'tt', 'frameRateHz');
%        
%    end
% end

%% 2P C1 sessions
% CC_paper_expts_crus;
% TP_img_dir = ['Y:\home\lindsey\Analysis\2P\Jake\'];
% bx_dir = ['Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\'];
% dest_folder_2P = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\two-photon C1\'];
% load([dest_folder, 'expNums_2P']);
% lobule_ID = 1;
% 
% for day_num = 1:3
%    if day_num ==1
%        curr_dataset = expt(day_num);
%        dest_id = 'day_1';
%    elseif day_num==2
%        curr_dataset = expt(day_num);
%        dest_id = 'day_N1';
%    elseif day_num ==3
%        curr_dataset = expt(day_num);
%        dest_id = 'day_N2';
%    end
%    sess_included =  expNums(day_num).mouse_name([~isnan(expNums(day_num).numDendrites(:,lobule_ID))],1); 
%    
%    for sess_num = 1:length(curr_dataset.mouse)
%        %select mouse
%        mouse = strtrim(curr_dataset.mouse(sess_num,:));
%        date = curr_dataset.date(sess_num,:);
%        this_session = [date, '_', mouse];
%        
%        %determine whether this mouse is included in this bx condition
%        mouse_rep = cell(size(sess_included,1),1);
%        mouse_rep(:) = {mouse};
%        mouse_ind = find(cellfun(@strcmp, sess_included, mouse_rep));
%        if isempty(mouse_ind)
%            continue
%        end
%        
%        %load mask and bx file
%        img_fn2 = [date '_img' mouse];
%        bx_data = get_bx_data(bx_dir, img_fn2);
%        load([TP_img_dir, this_session, '\', this_session, '_ROI_TCs'], 'mask3D');
%        load([TP_img_dir, this_session, '\', this_session, '_splitImage']);
%        load([TP_img_dir, this_session, '\', this_session, '_targetAlign.mat']);
%        
%        %Only keep the masks/TCs from Crus I
%        curr_lobs = curr_dataset.areas{sess_num};
%        if  strmatch(curr_lobs{1}, 'C1') & strmatch(curr_lobs{2}, 'C1')     %remove non-CrusI areas from targetAlign_events and targetAligndFoverF
%             crus_cells = find(maskCat ==1 | maskCat==2);%1 and 2
%        elseif strmatch(curr_lobs{1}, 'C1') & isempty(strmatch(curr_lobs{2}, 'C1'))
%            crus_cells = find(maskCat ==1);%1 and 2
%        elseif isempty(strmatch(curr_lobs{1}, 'C1')) & strmatch(curr_lobs{2}, 'C1')
%            crus_cells = find(maskCat==2);%1 and 2%2
%        end
%        targetAlign_events = targetAlign_events(:,[crus_cells],:);
%        targetAligndFoverF = targetAligndFoverF(:,[crus_cells],:);
%        
%        %copy and save data
%        save([dest_folder_2P, '\', dest_id, '\', this_session, '_masks'], 'mask3D', 'maskCat_map', 'split_img');
%        save([dest_folder_2P, '\', dest_id, '\', this_session, '_behavior'], 'bx_data');
%        save([dest_folder_2P, '\', dest_id, '\', this_session, '_time_courses'], 'targetAligndFoverF', 'maskCat', 'ind_omit', ...
%            'ind_omit_long', 'ind_omit_short', 'ind_rew', 'ind_rew_postomit', 'ind_rew_preomit', 'ind_unexp', 'ind_unexp_long', 'ind_unexp_short', ...
%            'postwin_frames', 'prewin_frames', 'tt', 'frameRateHz');
%        
%        %copy piezo data
%        if exist([TP_img_dir, this_session, '\', this_session, '_cueAlignPiezo.mat'])
%            current_source = [TP_img_dir, this_session, '\', this_session, '_cueAlignPiezo.mat'];
%            current_destination = [dest_folder_2P, dest_id, '\', this_session, '_piezo_data.mat'];
%            [copy_outcome, copy_msg] = copyfile(current_source, current_destination, 'f');
%            if copy_outcome == 0
%                list_of_failed_copies = [list_of_failed_copies; this_session];%save file name so it can be fixed later
%            end
%        end
%    end
% end

%% 2P C2 sessions

% CC_paper_expts_crus;
% TP_img_dir = ['Y:\home\lindsey\Analysis\2P\Jake\'];
% bx_dir = ['Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\'];
% dest_folder_2P = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\two-photon C2\'];
% load([dest_folder, 'expNums_2P']);
% lobule_ID = 2;
% 
% for day_num = 1:3
%    if day_num ==1
%        curr_dataset = expt(day_num);
%        dest_id = 'day_1';
%    elseif day_num==2
%        curr_dataset = expt(day_num);
%        dest_id = 'day_N1';
%    elseif day_num ==3
%        curr_dataset = expt(day_num);
%        dest_id = 'day_N2';
%    end
%    sess_included =  expNums(day_num).mouse_name([~isnan(expNums(day_num).numDendrites(:,lobule_ID))],1); 
%    
%    for sess_num = 1:length(curr_dataset.mouse)
%        %select mouse
%        mouse = strtrim(curr_dataset.mouse(sess_num,:));
%        date = curr_dataset.date(sess_num,:);
%        this_session = [date, '_', mouse];
%        
%        %determine whether this mouse is included in this bx condition
%        mouse_rep = cell(size(sess_included,1),1);
%        mouse_rep(:) = {mouse};
%        mouse_ind = find(cellfun(@strcmp, sess_included, mouse_rep));
%        if isempty(mouse_ind)
%            continue
%        end
%        
%        %load mask and bx file
%        img_fn2 = [date '_img' mouse];
%        bx_data = get_bx_data(bx_dir, img_fn2);
%        load([TP_img_dir, this_session, '\', this_session, '_ROI_TCs'], 'mask3D');
%        load([TP_img_dir, this_session, '\', this_session, '_splitImage']);
%        load([TP_img_dir, this_session, '\', this_session, '_targetAlign.mat']);
%        
%        %Only keep the masks/TCs from Crus II
%        curr_lobs = curr_dataset.areas{sess_num};
%        if  strmatch(curr_lobs{1}, 'C2') & strmatch(curr_lobs{2}, 'C2')     %remove non-CrusI areas from targetAlign_events and targetAligndFoverF
%             crus_cells = find(maskCat ==1 | maskCat==2);%1 and 2
%        elseif strmatch(curr_lobs{1}, 'C2') & isempty(strmatch(curr_lobs{2}, 'C2'))
%            crus_cells = find(maskCat ==1);%1 and 2
%        elseif isempty(strmatch(curr_lobs{1}, 'C2')) & strmatch(curr_lobs{2}, 'C2')
%            crus_cells = find(maskCat==2);%1 and 2%2
%        elseif isempty(strmatch(curr_lobs{1}, 'C2')) & isempty(strmatch(curr_lobs{2}, 'C2'))
%            continue
%        end
%        targetAlign_events = targetAlign_events(:,[crus_cells],:);
%        targetAligndFoverF = targetAligndFoverF(:,[crus_cells],:);
%        
%        %copy and save data
%        save([dest_folder_2P, dest_id, '\', this_session, '_masks'], 'mask3D', 'maskCat_map', 'split_img');
%        save([dest_folder_2P, dest_id, '\', this_session, '_behavior'], 'bx_data');
%        save([dest_folder_2P, dest_id, '\', this_session, '_time_courses'], 'targetAligndFoverF', 'maskCat', 'ind_omit', ...
%            'ind_omit_long', 'ind_omit_short', 'ind_rew', 'ind_rew_postomit', 'ind_rew_preomit', 'ind_unexp', 'ind_unexp_long', 'ind_unexp_short', ...
%            'postwin_frames', 'prewin_frames', 'tt', 'frameRateHz');
%        
%        %copy piezo data
%        if exist([TP_img_dir, this_session, '\', this_session, '_cueAlignPiezo.mat'])
%            current_source = [TP_img_dir, this_session, '\', this_session, '_cueAlignPiezo.mat'];
%            current_destination = [dest_folder_2P, dest_id, '\', this_session, '_piezo_data.mat'];
%            [copy_outcome, copy_msg] = copyfile(current_source, current_destination, 'f');
%            if copy_outcome == 0
%                list_of_failed_copies = [list_of_failed_copies; this_session];%save file name so it can be fixed later
%            end
%        end
%    end
% end

%% 2P LS piezo sessions

%'i098' 'i1032' 'i095' 'i096' 'i097'
CC_paper_expts_LS;   %need to replace with versions that only include the animals in the paper
TP_img_dir = ['Y:\home\lindsey\Analysis\2P\Jake\'];
bx_dir = ['Y:\home\jake\Data\WidefieldImaging\GCaMP\behavior\'];
dest_folder_2P = ['Y:\home\jake\Misc\classical conditioning paper\eLife paper data\two-photon LS piezo\'];
piezo_dir = ['Y:\home\jake\Data\2P_imaging\'];
load([dest_folder, 'expNums_LS_piezo']);

day_num = 2;
curr_dataset = expt(day_num);
dest_id = 'day_N1';
sess_included =  expNums.mouse_name([~isnan(expNums.numDendrites)],1); 

for sess_num = 16:length(curr_dataset.mouse);   %only include LS piezo mice
    %select mouse
    mouse = strtrim(curr_dataset.mouse(sess_num,:));
    date = curr_dataset.date(sess_num,:);
    this_session = [date, '_', mouse];
    
    %determine whether this mouse is included in this bx condition
    mouse_rep = cell(size(sess_included,1),1);
    mouse_rep(:) = {mouse};
    mouse_ind = find(cellfun(@strcmp, sess_included, mouse_rep));
    if isempty(mouse_ind)
        continue
    end
    
    %load mask and bx file
    img_fn2 = [date '_img' mouse];
    bx_data = get_bx_data(bx_dir, img_fn2);
    %load([TP_img_dir, this_session, '\', this_session, '_ROI_TCs'], 'mask3D');
    load([TP_img_dir, this_session, '\', this_session, '_targetAlign.mat']);
       
    %copy and save data
    %save([dest_folder_2P,  '\', dest_id, '\', this_session, '_mask3D'], 'mask3D');
    save([dest_folder_2P, '\',  dest_id, '\', this_session, '_behavior'], 'bx_data');
    save([dest_folder_2P, '\', dest_id, '\', this_session, '_time_courses'], 'targetAligndFoverF', 'ind_omit', ...
        'ind_omit_long', 'ind_omit_short', 'ind_rew', 'ind_rew_postomit', 'ind_rew_preomit', 'ind_unexp', 'ind_unexp_long', 'ind_unexp_short', ...
        'postwin_frames', 'prewin_frames', 'tt', 'frameRateHz');
     
    %copy piezo data
    current_source = [TP_img_dir, this_session, '\', this_session, '_cueAlignPiezo.mat'];
    current_destination = [dest_folder_2P, dest_id, '\', this_session, '_piezo_data.mat'];
    [copy_outcome, copy_msg] = copyfile(current_source, current_destination, 'f');
    if copy_outcome == 0
        list_of_failed_copies = [list_of_failed_copies; this_session];%save file name so it can be fixed later
    end
    
end

