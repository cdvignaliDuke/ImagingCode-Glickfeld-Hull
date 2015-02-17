function write_mov(img, params)

%params is the input parameters for writing the move:
% dest - movie destination
% speed - move speed up
%sampeling_rate - real sampling rate
% sz - dimension of each picture



aviobj = VideoWriter([params.dest '.avi']);

aviobj.FrameRate = params.sampeling_rate*params.speed; %b_data.input.frameImagingFrequencyHz;
open(aviobj);
cmap = colormap(params.colormap);

% ----- remove average
if(params.remove_avg || params.zscroe )
    avg_img = mean(img,2);    
    for i=1:size(img,2)
        img(:,i)  =  img(:,i) - avg_img;
    end    
end
% ----- divide by SD
if(params.zscroe)
    sd_img = std(img,[], 2);
    for i=1:size(img,2)
        img(:,i)  =  img(:,i)./sd_img;
    end    
end

% ------ set values between 0 and 1
max_val = max(img(:));
min_val = min(img(:));
img= (img - min_val)/(max_val - min_val);


 for i=1: size(img,2) 
      c_img = img(:,i);
      % ----- save as color movie       
      for k=1:size(c_img,1)
          rgb_inx =  ceil(c_img(k)*size(cmap,1));
          if(rgb_inx== 0)
              rgb_inx =1;
          end
          color_img(k,:) = cmap(rgb_inx,:);
      end
      
      save_img = [];
      for k=1:size(color_img,2)       
          
          re_shape = reshape(color_img(:,k), params.sz); 
          save_img(:,:,k) = re_shape;
      end
      writeVideo(aviobj,save_img);
      

 end
 close(aviobj);