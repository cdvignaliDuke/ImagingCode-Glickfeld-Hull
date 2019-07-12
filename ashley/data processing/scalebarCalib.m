function [x_um,y_um,sb_img_50um] = scalebarCalib(exptDateYYMMDD, obj, img, zm)%% for 30 Hz imaging on 2P
d = datestr(datenum(num2str(exptDateYYMMDD,'%d'),'yymmdd'),'yymmdd');
if strcmp(obj,'16x')
    if any(strcmp({'14','15','16'},d(1:2)))
        % 16x 2014 - 2016
        x_um = 555;
        y_um = 233;
    elseif any(strcmp({'17','18','19'},d(1:2)))
        % 16x 2017-2019
        x_um = 1030;
        y_um = 581;
    end
elseif strcmp(obj,'25x')
    if any(strcmp({'14','15'},d(1:2)))
        error('scale info not available for this date')
    elseif any(strcmp({'16','17','18','19'},d(1:2)))
        % 25x 2017-2019
        x_um = 673; 
        y_um = 382;
    end
end
if length(nargin > 3)
    x_um = x_um/zm;
    y_um = y_um/zm;
end


%%
if nargin > 2
    [y_pix,x_pix] = size(img);
    x_calib = x_um/x_pix;
    y_calib = y_um/y_pix;
    x50um = x_calib*50;
    y5um = y_calib*5;

    sb_x_ind = 500:500+x50um;
    sb_y_ind = 200:200+y5um;

    sb_img = zeros(y_pix,x_pix);
    sb_img(sb_y_ind,sb_x_ind) = 1;
    sb_img_50um = sb_img;
    figure;imagesc(sb_img)
end