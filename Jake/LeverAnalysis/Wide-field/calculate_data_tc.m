function calculate_data_tc(meta_data_dir, cluster_dir, bx_out_dir, stable_int);
%loads the raw data and uses that to generate a TC for each ROI. These TCs
%represent raw F and NOT a df/f. Saves data_tc and the frame size to the
%bx_outputs mat file. 

%load the meta_data and cluster info. meta_data used to load tiff files.
load(meta_data_dir, 'meta_data2');
load(cluster_dir);

%load tiff files
img = [];
for ii = 1:length(meta_data2)
    img_subset = readtiff(meta_data2{ii}(1).Filename);
    img = cat(3,img,img_subset);
    disp(['file #', num2str(ii), ' download complete']);
end
sz = size(img);
avg_img = mean(img(:,:,500:700),3);

%motion register the tiff
mean_img = squeeze(mean(mean(img,2),1))';
std_img = std(mean_img);
reg_ind = find(mean_img > std_img/2); 
img_ref = mean(img(:,:,[stable_int(1):stable_int(2)]),3);
[reg_out, img_reg] = stackRegister(img(:,:,[reg_ind]), img_ref);
img(:,:,[reg_ind]) = img_reg;
clear img_reg

%log the #of pixels in each ROI
roi_sz = sum(cluster.roi_mask,2);

% ----- cluster data by ROI and get raw TCs
disp('extracting ROIs separately');
data_tc = [];
for i = 1:cluster.num_cluster;
    data_tc(i,:) = stackGetTimeCourses(img, reshape(cluster.roi_mask(i,:),[sz(1) sz(2)]));
end

%save sz and data_tc to bx_outputs
save(bx_out_dir, 'data_tc', 'sz', 'avg_img', 'reg_out', 'reg_ind');
end




