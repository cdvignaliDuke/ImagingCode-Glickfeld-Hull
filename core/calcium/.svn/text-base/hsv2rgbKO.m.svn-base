function rgb = hsv2rgbKO (hueKO, saturation, lightness) 
%HSV2RGBKO 
% RGB = HSV2RGBKO (HUEKO, SATURATION, LIGHTNESS) 
% In hsv2rgb, red is opposed to cyan.
% In hsv2rgbKO, red is opposed to green, and blue is opposed to yellow
% Kenichi Ohki 2004.9.23
% Kenichi Ohki 2008.11.5 updated to fix a problem of lower lightness at purple and cyan.
%                           also now this function can accept any matrix.

rgb=hsv2rgb(hueKO2HSV(hueKO),saturation,lightness);


% white = [1;1;1];
% red = [1;0;0];
% yellow = [1;1;0];
% green= [0;1;0];
% blue = [0;0;1];
% 
% if hue < 0.25
%     rgb = red * (1-hue*4) + yellow * hue*4;
% else
%     if hue < 0.5
%         hue = hue - 0.25;
%         rgb = yellow * (1-hue*4) + green * hue*4;
%     else if hue < 0.75
%             hue = hue -0.5;
%             rgb = green * (1-hue*4) + blue * hue*4;
%         else
%             hue = hue -0.75;
%             rgb = blue * (1-hue*4) + red * hue*4;
%         end
%     end
% end
% rgb= rgb.*saturation + white.*(1-saturation);
% rgb= rgb.*lightness;
% 
% 
                