clear
days = {'151021_img29', '151022_img29', '151009_img30', '151011_img30', '151211_img32', '151212_img32', '160129_img35', '160131_img35', '160129_img36', '160131_img36', '160314_img38', '160315_img38', '160319_img41', '160320_img41', '160606_img46', '160722_img53', '160904_img55'}; %'150718_img27', '150719_img27', '150716_img28', '150717_img28', 
bx_source          = ['Z:\Data\WidefieldImaging\GCaMP\behavior\'];
image_source_base  = ['Z:\Data\WidefieldImaging\GCaMP\']; %location of permanently stored image files for retreiving meta data
image_dest_base    = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\']; %stores the data on crash in the lever analysis folder
bx_outputs_dir     = ['Z:\Analysis\WF Lever Analysis\BxAndAnalysisOutputs\BxOutputs\'];
PCA_output_dir_base= ['Z:\Analysis\WF Lever Analysis\PCA_output_dir\'];
kmeans_output_dir_base  = ['Z:\Analysis\WF Lever Analysis\kmeans_output_dir\'];
old_cd = cd; %save old cd so I can restore it later
downsample_factor = 0.2;

for ii = [15, 16, 17]
    %pauses inserted to allow for the viewing and selection of PCs
    days{ii}
    PCA_output_dir = [PCA_output_dir_base days{ii} '\'];
    kmeans_output_dir = [kmeans_output_dir_base days{ii} '\'];

%     combined_release_frames = readtiff([PCA_output_dir, days{ii}, '_combined_release_movie.tif']);
%     combined_release_frames = double(imresize(combined_release_frames, downsample_factor));
%     save([PCA_output_dir 'down sampled cropped epoch specific movies\', days{ii}, '_combined_release_movie'], 'combined_release_frames');
%     
%     lick_frames = readtiff([PCA_output_dir, days{ii}, '_licking_movie.tif']);
%     lick_frames = double(imresize(lick_frames, downsample_factor));
%     save([PCA_output_dir 'down sampled cropped epoch specific movies\', days{ii}, '_lick_frames'], 'lick_frames');
%     
%     iti_frames = readtiff([PCA_output_dir, days{ii}, '_iti_movie.tif']);
%     iti_frames = double(imresize(iti_frames, downsample_factor));
%     save([PCA_output_dir 'down sampled cropped epoch specific movies\', days{ii}, '_iti_frames'], 'iti_frames');
    
%     early_release_frames = readtiff([PCA_output_dir, days{ii}, '_early_release_movie.tif']);
%     early_release_frames = double(imresize(early_release_frames, downsample_factor));
%
%     corr_release_frames = readtiff([PCA_output_dir, days{ii}, '_corr_release_movie.tif']);

    load([PCA_output_dir, days{ii}, '_corr_full_movie']);
    corr_full_frames = double(imresize(corr_full_frames, downsample_factor));
    save([kmeans_output_dir , days{ii}, '_corr_full_movie'], 'corr_full_frames');

    load([PCA_output_dir, days{ii}, '_early_full_movie']);
    early_full_frames = double(imresize(early_full_frames, downsample_factor));
    save([kmeans_output_dir, days{ii}, '_early_full_movie'], 'early_full_frames');


end