function [ROI_x ROI_y] =  get_movie_ROI(meta_data2, frame_info)

image_dest = meta_data2{1}(1).Filename;
ROI_x  =[];
ROI_y = [];

% read 10 frames across the movie
NUM_FRAMES = 20;
first_frame = 1; 
last_frame = length(meta_data2{1});
use_inx = ceil(linspace(first_frame, last_frame, NUM_FRAMES));
for i=1:length(use_inx)
    c_img =  double(imread(image_dest, use_inx(i)));
    img(:,:,i) = c_img;
end

% calculate STD and mean
sd = std(img,[],3);
avg = mean(img,3);

done =0;
roi_y = [1 size(sd,1)]; %sd is the size of one full frame
roi_x =  [1 size(sd,2)];
figure;
while(~done)
    % plot waveforms
    hold off;
    pcolor(avg); axis ij; shading flat; hold on; colormap jet;
    % plot square
    plot( [roi_x(1) roi_x(1)], roi_y, 'w', 'linewidth', 2);
    plot( [roi_x(2) roi_x(2)], roi_y, 'w', 'linewidth', 2);
    plot( roi_x, [roi_y(1), roi_y(1)], 'w', 'linewidth', 2);
    plot( roi_x, [roi_y(2), roi_y(2)], 'w', 'linewidth', 2);
     
    % get top right corner
    axis tight;
    title('Select top left then bottom right corner. Press s when completed', 'FontSize' ,20);
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
    title('Select top left then bottom right corner. Press s when completed', 'FontSize' ,20);
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


