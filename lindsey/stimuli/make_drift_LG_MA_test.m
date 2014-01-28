%make a square wave grating, save it.. 


LUM_HI = 0; %.95; %set to 0 and 1 for prev version
LUM_LO = 1; %.4;
LUM_MEAN = .7565; %figure this out impirically using calibration
                  %using full field blanks (see end of this file,
                  %after the word 'break', 
                  %for how to make these)
		  
scale1 = 8; %scales the pixelation so that movies are smaller and
            %run better. this is good to do for square wave
            %movies.. 
	    
%convert degrees into pixels: 
%at 20cm, 1cycle = 20 degrees => tan(20)=O/20
%so O = 20*tan(20*pi/180) %in cm
pix_per_cm = 24./scale1;

%Ly = 1024/scale1;
%Lx = 768/scale1;
Ly = 1280/scale1; %choose appropriate values for your monitor
Lx = 1024/scale1;

%ori_vec = [0 45 90 135]*pi/180;
%ori_vec = [0 90]*pi/180; %[0 90]*pi/180;
%ori_vec = [90 0 75 60 45 30 15]*pi/180; %[0 90]*pi/180;


%ori_vec = [90 0 10 22.5 45 67.5 80]*pi/180; %[0 90]*pi/180;
%ori_vec = 0*pi/180; %[0 45 90 135 180 225 270 315]*pi/180; % 0 26 45 58 67.5 74 79 82 84.5 86]*pi/180; %[0 90]*pi/180;

%ori_vec = [0 45 90 135 180 225 270 315]*pi/180;
ori_vec = [0 45]*pi/180;  

PWD_visfiles = ['\\Zmey\storlab\users\Lindsey\visual stimuli\'];
	

a = dir(PWD_visfiles);
if isempty(a)
  eval(['!mkdir ',PWD_visfiles]);
end

Dist = 18; %cm from screen 
f = .1; %cycles per degree
Contrast = 100; %percent
speed = 15; %Hz %DriftRate*Nphase;
DriftRateMat =1; %[0 1 2 3]; %#periods/sec
DriftRate = DriftRateMat(1);
Nphase = speed./DriftRate; %makwe sure tihs is a whole number
    

Ncyc = 2;
MAKE_BLANK_POST = 1; %1 -> make same length movie of blank frames added to end.
 
color = 'w'; %can be 'b' or 'w' right now.. 
	     %now make a 3color version: 
if strcmp(color,'w')
  Cmap1 = [1 1 1];
elseif strcmp(color,'b');
  Cmap1 = [0 0 1];
end
 
%for now, assume screen is 480x640


phase = 0*pi/180; %phase of grating (in rad)

X = [1:Lx]'*ones(1,Ly);
Y = ones(Lx,1)*[1:Ly];

Opp =  Dist*tan((1/f)*pi/180)* pix_per_cm %length of 1 cycle, in pixels

clear('M_ALL');
step_ALL = 0;
for count_ori = 1:length(ori_vec)
  ori = ori_vec(count_ori);
  ori_txt = num2str(round(ori*180/pi*100));
  if length(ori_txt) == 1
    ori_txt = ['000',ori_txt];
  elseif length(ori_txt) == 2
    ori_txt = ['00',ori_txt];
  elseif length(ori_txt) == 3
    ori_txt = ['0',ori_txt];
  end
  
  oriB = round(ori*180/pi);
  %first, generate a sine wave grating..
  Big_Z3 = zeros(Nphase,Lx,Ly,3);
  for count_phase = 1:Nphase

      phase = count_phase./Nphase*2*pi;
      %elseif DriftRate = 0;
      %  phase = 0;
      %end

      Z = sin((X*sin(ori) + Y*cos(ori))./Opp + phase);
      Z2a = Z>0;
      %  imagesc(Z2); colormap('gray'); axis image; pause(.2)


      %Z2_use = (Z2a-.5)*Contrast/100 + .5;
      %assume Z2a is all 0's and 1's
      Z2_use = zeros(size(Z2a));
      Z2_use(find(Z2a==1)) = LUM_HI;
      Z2_use(find(Z2a==0)) = LUM_LO;


      Z3 = zeros(Lx,Ly,3);
      for count = 1:3
          Z3(:,:,count) = Z2_use*Cmap1(count);
      end
      image(Z3); axis image; pause(.2)

      Big_Z3(count_phase,:,:,:) = Z3;
  end

  %now make an avi
  Nframes = Ncyc * Nphase;
  avi_Z3 = zeros(Nframes,Lx,Ly,3);

 
  %for count_drift = 1:length(DriftRateMat)
      clear M
      count_phase2a = 0;
      for n0 = 1:Ncyc
          for count_phase = 1:Nphase
              count_phase2a = count_phase2a + 1;

              if DriftRate == 1
                  count_phase2 = rem(count_phase2a,Nphase) + 1;
              elseif DriftRate == 0
                  count_phase2 = 1;
              elseif DriftRate == 2
                  %slow down by a factor of 2:
                  count_phase2 = rem(ceil(count_phase2a/2),Nphase) + 1
              elseif DriftRate == 3
                  count_phase2 = rem(ceil(count_phase2a/3),Nphase) + 1
                  %	  count_phase2 = count_phase;
              end


              Z3B = squeeze(Big_Z3(count_phase2,:,:,:));
              h = figure(gcf); clf;

              h2 = image(Z3B); axis image;
              axis off;
              pause(.1);
              %M((n0-1)*2 + count_rev) = getframe(gcf,[440 314 560 420]);;
              M((n0-1)*Nphase + count_phase) = im2frame(Z3B);
              step_ALL = step_ALL + 1;
              M_ALL(step_ALL) = im2frame(Z3B);
          end
      end
      Nframes_tot = (n0-1)*Nphase + count_phase;
      if MAKE_BLANK_POST == 1

          Z4 = zeros(Lx,Ly,3);
          for count = 1:3
              Z4(:,:,count) = LUM_MEAN*Cmap1(count);
          end

          count_phase2a = 0;
          for n0 = 1:Ncyc
              for count_phase = 1:Nphase
                  count_phase2a = count_phase2a + 1;

                  if DriftRate == 1
                      count_phase2 = rem(count_phase2a,Nphase) + 1;
                  elseif DriftRate == 0
                      count_phase2 = 1;
                  elseif DriftRate == 2
                      %slow down by a factor of 2:
                      count_phase2 = rem(ceil(count_phase2a/2),Nphase) + 1
                  elseif DriftRate == 3
                      count_phase2 = rem(ceil(count_phase2a/3),Nphase) + 1
                      %	  count_phase2 = count_phase;
                  end


                  %Z3B = squeeze(Big_Z3(count_phase2,:,:,:));
                  Z3B = Z4;

                  h = figure(gcf); clf;

                  h2 = image(Z3B); axis image;
                  axis off;
                  pause(.1);
                  %M((n0-1)*2 + count_rev) = getframe(gcf,[440 314 560 420]);;
                  M(Nframes_tot + (n0-1)*Nphase + count_phase) = im2frame(Z3B);
                  step_ALL = step_ALL + 1;
                  M_ALL(step_ALL) = im2frame(Z3B);
              end
          end

      end

  file_save0 = [PWD_visfiles,'090804_vid_drift_grat_4dir_', ...
	       num2str(Nphase),'_or1_',ori_txt,'_c',color, ...
	       '_P',num2str(1/f),'_cont_',num2str(Contrast),'B', ...
		'_DR_',num2str(DriftRate),'_Lxy_',num2str(Lx),'_',num2str(Ly),'_LUM_',num2str(LUM_LO),'_',num2str(LUM_HI),'.avi'];
  movie2avi(M,file_save0,'compression','none','fps',speed);
%  movie2avi(M3,file_save0,'compression','none');
  %end
  
end

ori_txt_all = [];
    
for count = 1:length(ori_vec)
    ori_txt_all = [ori_txt_all,'_',num2str(ori_vec(count)*180/pi)];
end

  file_save0 = [PWD_visfiles,'090804_ALLORI_vid_drift_grat_4dir_', ...
	       num2str(Nphase),'_or1_',ori_txt_all,'_c',color, ...
	       '_P',num2str(1/f),'_cont_',num2str(Contrast),'B', ...
		'_DR_',num2str(DriftRate),'_Lxy_',num2str(Lx),'_',num2str(Ly),'_LUM_',num2str(LUM_LO),'_',num2str(LUM_HI),'.avi'];
  movie2avi(M_ALL,file_save0,'compression','none','fps',speed);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


break

%make a blank screen with the medium color; 

mean_lum = .7565; %.7565 from June20.7; %.35;
Z4 = zeros(Lx,Ly,3);
for count = 1:3
  Z4(:,:,count) = mean_lum*Cmap1(count);
end

strP = ['_Ppt',num2str(mean_lum*100),'_'];


file_save = [PWD_visfiles,'Jmeanlum',strP,'_c',color,'.jpg'];
imwrite(Z4,file_save,'bitdepth',8,'Quality',100);

file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'.tif'];
imwrite(Z4,file_save,'compression','none');

%make an avi of teh blank screen
file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'.avi'];
M =  im2frame(Z4);
movie2avi(M,file_save,'compression','None');



%make files of many luminances: 

for mean_lum = .3:.01:1
  Z4 = zeros(Lx,Ly,3);
  for count = 1:3
    Z4(:,:,count) = mean_lum*Cmap1(count);
  end
  
  strP = ['_Ppt',num2str(mean_lum*100),'_'];
  
  file_save = [PWD_visfiles,'Jmeanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.jpg'];
  imwrite(Z4,file_save,'bitdepth',8,'Quality',100);

  file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.tif'];
imwrite(Z4,file_save,'compression','none');
  
  
  %make an avi of teh blank screen
  file_save = [PWD_visfiles,'meanlum',strP,'_c',color,'_Lxy_',num2str(Lx),'_',num2str(Ly),'.avi'];
  M =  im2frame(Z4);
  movie2avi(M,file_save,'compression','None');
  
end

