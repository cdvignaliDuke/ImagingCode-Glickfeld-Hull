%% for 30 Hz imaging on 2P
% 16x 2014 - 2016
x_um = 555;
y_um = 233;

% 25x 2017
x_um = 673; 
y_um = 382;

% 16x 2018
x_um = 1030;
y_um = 581;

%%
[y_pix,x_pix] = size(img);
x_calib = x_um/x_pix;
y_calib = y_um/y_pix;
x50um = x_calib*50;
y5um = y_calib*5;

sb_x_ind = 500:500+x50um;
sb_y_ind = 200:200+y5um;

sb_img = zeros(y_pix,x_pix);
sb_img(sb_y_ind,sb_x_ind) = 1;

figure;imagesc(sb_img)
