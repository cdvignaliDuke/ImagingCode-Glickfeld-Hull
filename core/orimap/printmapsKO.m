function p=printmapsKO(print_flag,inpInt,inpStr)
% modeled from domaplistSY.m
p=0;
debugflag=false;

if nargin < 1
    print_flag=2;
end
if nargin<2
    inpInt=0;
end
if nargin<3
    inpStr='N/A';
end

s=what;
dirname=s.path;



h=figure;
set(h,'PaperUnits','normalized');
set(h,'PaperPosition',[0.05 0.05 0.9 0.9]);

%main comment
subplot('Position',[0.05 0.95 0.9 0.05]); 
set(gca,'FontSize',6);
axis off;
text(0.1,0.8,[num2str(inpInt)],'FontSize',6);
text(0.1,0.6,inpStr,'FontSize',6);
text(0.1,0.4,dirname,'FontSize',6);



% images
%manualy place each image on page
imtype='FOV.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,9);
    axes('units','normal','Position', [0.7 0.5 0.09 0.09]);
    imshow(im);
    set(gca,'FontSize',6);
    title(strrep(imtype,'_','\_'));
    
    for i=0:7
    axes('units','normal','Position', [0.2+0.1*i 0.85 0.09 0.09]);
    imshow(im);
end
    
    for i=0:7
    axes('units','normal','Position', [0.1 0.75-0.1*i 0.09 0.09]);
    imshow(im);
end
    
    
    
    
end


imtype='dir_angle.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
   % subplot(11,8,10);
    axes('units','normal','Position', [0.8 0.5 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_angle.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,11);
        axes('units','normal','Position', [0.9 0.5 0.09 0.09]);

    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end



for i=1:8
    imtype=['dir_dF_' num2str(i) '.tif'];
    imname=[dirname '\' imtype];
    if exist(imname)>0
        im=imread(imname);
   %     subplot(11,8,24+i);
        axes('units','normal','Position', [0.2+(i-1)*0.1 0.75 0.09 0.09]);
        set(gca,'FontSize',6);
        imshow(im);
    title(strrep(imtype,'_','\_'));
    end
end

for i=1:8
    imtype=['dir_ratio_' num2str(i) '.tif'];
    imname=[dirname '\' imtype];
    if exist(imname)>0
        im=imread(imname);
  %      subplot(11,8,32+i);
        axes('units','normal','Position', [0.2+(i-1)*0.1 0.65 0.09 0.09]);
        set(gca,'FontSize',6);
        imshow(im);
    title(strrep(imtype,'_','\_'));
    end
end

for i=1:4
    imtype=['ori_dF_' num2str(i) '.tif'];
    imname=[dirname '\' imtype];
    if exist(imname)>0
        im=imread(imname);
 %       subplot(11,8,56+i);
         axes('units','normal','Position', [0.2+(i-1)*0.1 0.55 0.09 0.09]);

        set(gca,'FontSize',6);
        imshow(im);
    title(strrep(imtype,'_','\_'));
    end
end

for i=1:4
    imtype=['dir_ratio_' num2str(i) '.tif'];
    imname=[dirname '\' imtype];
    if exist(imname)>0
        im=imread(imname);
 %       subplot(11,8,64+i);
        axes('units','normal','Position', [0.2+(i-1)*0.1 0.45 0.09 0.09]);
        set(gca,'FontSize',6);
        imshow(im);
    title(strrep(imtype,'_','\_'));
    end
end


imtype='dir_dF_HLS.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,18);
        axes('units','normal','Position', [0.2 0.35 0.09 0.09]);

    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='dir_dF_HLS_hc.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,19);
        axes('units','normal','Position', [0.3 0.35 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='dir_dF_polar.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
  %  subplot(11,8,22);
        axes('units','normal','Position', [0.4 0.35 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='dir_dF_polar_hc.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,23);
            axes('units','normal','Position', [0.5 0.35 0.09 0.09]);

    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='dir_ratio_HLS.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
 %   subplot(11,8,42);
            axes('units','normal','Position', [0.6 0.35 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='dir_ratio_HLS_hc.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
 %   subplot(11,8,43);
            axes('units','normal','Position', [0.7 0.35 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='dir_ratio_polar.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,46);
            axes('units','normal','Position', [0.8 0.35 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='dir_ratio_polar_hc.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
 %   subplot(11,8,47);
            axes('units','normal','Position', [0.9 0.35 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end







imtype='ori_dF_HLS.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,18);
        axes('units','normal','Position', [0.2 0.25 0.09 0.09]);

    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_dF_HLS_hc.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,19);
        axes('units','normal','Position', [0.3 0.25 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='ori_dF_polar.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
  %  subplot(11,8,22);
        axes('units','normal','Position', [0.4 0.25 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_dF_polar_hc.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,23);
            axes('units','normal','Position', [0.5 0.25 0.09 0.09]);

    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_ratio_HLS.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
 %   subplot(11,8,42);
            axes('units','normal','Position', [0.6 0.25 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_ratio_HLS_hc.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
 %   subplot(11,8,43);
            axes('units','normal','Position', [0.7 0.25 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_ratio_polar.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
%    subplot(11,8,46);
            axes('units','normal','Position', [0.8 0.25 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_ratio_polar_hc.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
 %   subplot(11,8,47);
            axes('units','normal','Position', [0.9 0.25 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end




imtype='dir_dF_ave_change.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
 %   subplot(11,8,17);
            axes('units','normal','Position', [0.2 0.15 0.09 0.09]);

    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='dir_dF_mag.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.3 0.15 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='dir_dF_max_change.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.4 0.15 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='dir_dF_tune.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.5 0.15 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end




imtype='dir_ratio_ave_change.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.6 0.15 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='dir_ratio_mag.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.7 0.15 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='dir_ratio_max_change.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.8 0.15 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='dir_ratio_tune.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.9 0.15 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end



imtype='ori_dF_ave_change.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
 %   subplot(11,8,17);
            axes('units','normal','Position', [0.2 0.05 0.09 0.09]);

    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='ori_dF_mag.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.3 0.05 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_dF_max_change.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.4 0.05 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='ori_dF_tune.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.5 0.05 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end




imtype='ori_ratio_ave_change.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.6 0.05 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='ori_ratio_mag.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.7 0.05 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end

imtype='ori_ratio_max_change.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.8 0.05 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end


imtype='ori_ratio_tune.tif';
imname=[dirname '\' imtype];
if exist(imname)>0
    im=imread(imname);
            axes('units','normal','Position', [0.9 0.05 0.09 0.09]);
    set(gca,'FontSize',6);
    imshow(im);
    title(strrep(imtype,'_','\_'));
end








if print_flag==1
    print;
end

% output as a postscript file
% added by KO 09/15/04

if print_flag==2
    print  %YC 09/17/04
    print ('-dpsc','-append',[dirname,'_map_plot.ps']);
end


warning on all   

% if debugflag
%     disp('Pause,  hit CR');
%     pause
% end
% close all

