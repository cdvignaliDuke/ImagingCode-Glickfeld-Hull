ICbad = [];

kl = floor(nIC/5);

colord=[         0         0    1.0000
    0    0.4000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.8, 0.5, 0
    0         0    0.5
    0         0.85      0];

for k = 1:kl
    fig = figure;
    shift  = 0;
    for i = (k-1)*5+1:k*5
        plot(tc_avg(:,i)+shift,'Color',colord(mod(i-1,size(colord,1))+1,:)); hold on
        %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
        %         hline(shift);
        %     end
        shift = shift+10000;
    end
    set(gca,'YTick',10000:10000:shift);
    set(gca,'YTicklabel',((k-1)*5+1:k*5));
end

if mod(nIC,5) ~= 0
    figure;
    shift  = 0;
    for i = 5*kl+1:nIC
        plot(tc_avg(:,i)+shift,'Color',colord(mod(i-1,size(colord,1))+1,:)); hold on
        %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
        %         hline(shift);
        %     end
        shift = shift+10000;
    end
    set(gca,'YTick',10000:10000:shift);
    set(gca,'YTicklabel',(5*kl+1:nIC));
end

shift  = 0;
    fig = figure;
    for i = 1:10
        plot(npSubTC(:,i)+shift,'Color',colord(mod(i-1,size(colord,1))+1,:)); hold on
        %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
        %         hline(shift);
        %     end
        shift = shift+10000;
    end
    set(gca,'YTick',10000:10000:shift);
    set(gca,'YTicklabel',(1:nCells));
   
file_info;
sub = 9;
rID = 1;
mouse = mouseID{sub};
dateID = date;
date = dateID{sub};
dest_sub=fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[dateID{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
cd(dest_sub)
load('Results.mat');
load('ROI_TCs.mat');
clear input
ICbad_input = input('Number of bad IC ', 's');
ICbad = [str2num(ICbad_input)];
tc_avg(:,ICbad) = []; mask3D(:,:,ICbad) = [];


save( 'Results.mat', 'tc_avg', 'mask_raw', 'mask_flat', 'mask3D', 'mask_final', 'data_corr');
save( 'ROI_TCs.mat','tc_avg', 'mask_final', 'sz');
close all