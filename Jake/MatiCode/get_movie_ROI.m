function [ROI_x ROI_y] =  get_movie_ROI(moive_info, image_dest, first_frame, last_frame)
ROI_x  =[];
ROI_y = [];
% read 10 frames across the movie
NUM_FRAMES = 10;
first = first_frame;
last =last_frame;

use_inx = ceil(linspace(first, last, NUM_FRAMES));

for i=1:length(use_inx)
    
      c_img =  double(imread(image_dest, use_inx(i), 'info', moive_info));
     img(:,:,i) = c_img;
end
% calculate STD and plot on figure;
sd = std(img,[],3);
avg = mean(img,3);
%figure; pcolor(sd); axis ij; shading flat; %colormap gray;

done =0;
roi_y = [1 size(sd,1)];
roi_x =  [1 size(sd,2)];
% get ROI  press 'd'
figure; 
while(~done)    
    % plot waveforms
    hold off;
     %pcolor(sd); axis ij; shading flat; hold on;%colormap gray;
     pcolor(avg); axis ij; shading flat; hold on;%colormap gray;
    % plot square
    plot( [roi_x(1) roi_x(1)], roi_y, 'w', 'linewidth', 2);
    plot( [roi_x(2) roi_x(2)], roi_y, 'w', 'linewidth', 2);
    plot( roi_x, [roi_y(1), roi_y(1)], 'w', 'linewidth', 2);
    plot( roi_x, [roi_y(2), roi_y(2)], 'w', 'linewidth', 2);
    
    % get top right corner
    axis tight;
    title('press Top LEFT corner or S to stop', 'FontSize' ,20);
    [ x_tr , y_tr, butt] = ginput(1);
    if(butt ==1)
        
        if(x_tr > size(sd,1))
            x_tr= size(sd,1);
        end
        if(y_tr > size(sd,2))
            y_tr= size(sd,2);
        end
    elseif(char(butt) == 's')
         done =1;
         break;
    end
    
    % get buttom left corner
    title('press Bottom RIGHT corner or s to stop', 'FontSize' ,20);
    [ x_bl , y_bl, butt] = ginput(1);
    if(butt ==1)
        
        if(x_bl<1)
            x_bl= 1;
        end
        if(y_bl <1)
            y_bl= 1;
        end
    elseif(char(butt) == 's')
        done =1;
        break;
    end
    roi_x = round([ x_bl x_tr]);
    roi_y = round([y_bl y_tr]);
    
end
% switch x and y - because in pcolor x and y are switched
if(roi_y(1) > roi_y(2))
    ROI_x = roi_y(2):roi_y(1);   
else
    ROI_x = roi_y(1):roi_y(2);   
end
    
if(roi_x(1) > roi_x(2))
    ROI_y = roi_x(2):roi_x(1);   
else
    ROI_y = roi_x(1):roi_x(2);   
end
%ROI_y = roi_x(1):roi_x(2);

    
