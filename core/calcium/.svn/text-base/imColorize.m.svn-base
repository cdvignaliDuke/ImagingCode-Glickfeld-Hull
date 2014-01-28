function out = imColorize (image, color)

% colorize gray image with a color
% image should be between 0-1
% color also should be between 0-1
% 2008. 6. 19. Kenichi Ohki 

out(:,:,1)=image.*color(1);
out(:,:,2)=image.*color(2);
out(:,:,3)=image.*color(3);

out(find(out<0))=0;
out(find(out>1))=1;
