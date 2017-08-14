

clear;
days = {'170321_img86', '170322_img86', '170323_img86', '170325_img86', '170326_img86', '170327_img86', '170330_img86'};
WF_CRP_licking_only_overview;

clear;
days = {'170408_img87', '170410_img87', '170411_img87', '170412_img87', '170413_img87', '170414_img87', '170415_img87', '170417_img87', '170418_img87', '170419_img87', '170420_img87', '170422_img87', '170423_img87', '170424_img87', '170425_img87', '170426_img87'};
WF_CRP_licking_only_overview;

clear;
days = {'170408_img88', '170410_img88', '170411_img88', '170412_img88', '170413_img88', '170414_img88', '170415_img88', '170417_img88', '170418_img88', '170419_img88', '170420_img88', '170421_img88', '170422_img88', '170423_img88', '170424_img88', '170425_img88', '170426_img88'};
WF_CRP_licking_only_overview;

clear;
days = {'170416_img90', '170417_img90', '170418_img90', '170419_img90', '170422_img90', '170423_img90', '170424_img90', '170425_img90', '170426_img90', '170427_img90', '170428_img90', '170429_img90', '170501_img90', '170502_img90', '170503_img90', '170504_img90', '170505_img90', '170506_img90', '170508_img90', '170509_img90', '170510_img90', '170511_img90', '170512_img90', '170513_img90', '170515_img90'};
WF_CRP_licking_only_overview;

clear;
days = {'170417_img91', '170418_img91', '170419_img91', '170420_img91', '170422_img91', '170423_img91', '170424_img91', '170425_img91', '170426_img91', '170427_img91', '170428_img91', '170429_img91', '170501_img91', '170502_img91', '170503_img91', '170504_img91'};
WF_CRP_licking_only_overview;

clear;
days = {'170420_img92', '170422_img92', '170423_img92', '170424_img92', '170425_img92', '170426_img92', '170427_img92', '170428_img92', '170429_img92', '170501_img92', '170503_img92', '170504_img92', '170505_img92', '170506_img92', '170509_img92'};
WF_CRP_licking_only_overview;

clear;
days = {'170510_img93', '170511_img93', '170512_img93', '170513_img93', '170515_img93', '170516_img93', '170517_img93', '170518_img93', '170519_img93', '170520_img93', '170522_img93', '170523_img93', '170524_img93', '170525_img93'};
WF_CRP_licking_only_overview;

clear;
days = {'170513_img89', '170515_img89', '170516_img89', '170517_img89', '170518_img89', '170519_img89', '170522_img89', '170523_img89', '170524_img89', '170525_img89', '170526_img89', '170527_img89', '170529_img89', '170530_img89', '170531_img89'};
WF_CRP_licking_only_overview;

clear;
days = {'170524_img94', '170525_img94', '170526_img94', '170527_img94', '170529_img94', '170530_img94', '170531_img94', '170601_img94', '170602_img94', '170604_img94', '170605_img94', '170606_img94'};
WF_CRP_licking_only_overview;

clear;
days = {'170605_img95', '170606_img95', '170607_img95', '170608_img95', '170609_img95', '170610_img95', '170611_img95', '170612_img95', '170613_img95', '170614_img95'};
WF_CRP_licking_only_overview;

clear;
days = {'170612_img96', '170613_img96', '170614_img96', '170615_img96', '170620_img96', '170621_img96', '170622_img96', '170623_img96', '170624_img96'};
WF_CRP_licking_only_overview;

clear;
days = {'170628_img98', '170629_img98', '170701_img98', '170703_img98', '170704_img98', '170705_img98', '170706_img98'};
WF_CRP_licking_only_overview;

clear;
days = {'170705_img99', '170706_img99', '170707_img99', '170708_img99', '170710_img99', '170711_img99', '170712_img99'};
WF_CRP_licking_only_overview;













%%
% clear
% tic
% day = '151009_img30';
% %day = '160606_img46';
% image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
% image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
% image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
% bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
% if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
%     mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
% end
% old_cd = cd; %save old cd so I can restore it later
% [all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
% frame_times = get_frame_time_by_movie_info(meta_data);
% dest =  [image_dest '\' day '_ROI'];
% save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
% save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
% shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
% toc
% clear
% %%
% clear
% tic
% day = '151011_img30';
% %day = '160606_img46';
% image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
% image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
% image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
% bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
% if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
%     mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
% end
% old_cd = cd; %save old cd so I can restore it later
% [all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
% frame_times = get_frame_time_by_movie_info(meta_data);
% dest =  [image_dest '\' day '_ROI'];
% save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
% save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
% shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
% toc
% clear
% %%
% clear
% tic
% day = '151021_img29';
% %day = '160606_img46';
% image_source = ['Z:\Data\WidefieldImaging\GCaMP\', day]; %location of permanently stored image files for retreiving meta data
% image_source_temp = ['C:\Users\jake\TempData\', day]; %looks at temporary storage to find the raw image data
% image_dest   = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\' day]; %stores the data on crash in the lever analysis folder
% bx_source    =  'Z:\Data\WidefieldImaging\GCaMP\behavior\';
% if exist(image_dest, 'file') ~= 7;   %check to make sure that image_dest exists as a folder in the correct location. 
%     mkdir(image_dest);               %if it does not exist then create that folder. Should make this into a function and insert into the 2P code. 
% end
% old_cd = cd; %save old cd so I can restore it later
% [all_files, meta_data, meta_data2] = obtain_tif_meta_data(image_source);
% frame_times = get_frame_time_by_movie_info(meta_data);
% dest =  [image_dest '\' day '_ROI'];
% save([dest '_frame_times'],  'frame_times');  %ALTER PATH so it saves in WF Lever Analysis
% save([dest '_meta_data'],  'meta_data', 'meta_data2', 'all_files');
% shrink_movie(day, image_source_temp, image_dest, meta_data2); %automatically saves the shrunken movie
% toc
% clear

