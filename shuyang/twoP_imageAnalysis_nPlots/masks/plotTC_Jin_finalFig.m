%used in twoP_TC_n_speed,this function is the same as plotTC_Jin, but
%better appearance and will save the plot to the final figure folder.

function plotTC_Jin_finalFig(x,tc_avg, sm, reg_mean, ICuse, frGetFrameRate)
%nmask = max(max(mask_final));
% RGBmat = NaN(nmask,3);
%
% c = [.1, .5, .9]; c2 = NaN(1, nmask);
% for cc = 1:3:nmask
%     c2(1,[cc:(cc+2)]) = c;
% end
%
%
%
% for a = 1:size(RGBmat,1);
%     one = a/nmask;
%     two = c2(1,a);
%     three = (nmask-a)/nmask;
%     RGBmat(a,:) = [one, two, three];
% end
% timeCourses = tc_avg;
% delta = double(5*nanmean(nanstd(timeCourses)));
% [nsamples, ncells] = size(timeCourses);
% offsets = repmat(0:delta:(ncells-1)*delta,nsamples,1);

% tt = [0:size(tc_avg,1)-1]/frGetFrameRate;
% figure;
% for ab = 1:size(RGBmat,1);
%     hold on
%     fig = plot(tt,timeCourses(:,ab)+offsets(:,ab),'color',RGBmat(ab,:));
% end
colord=[         0         0    1.0000
    0    0.4000         0
    1.0000         0         0
    0    0.7500    0.7500
    0.7500         0    0.7500
    0.8, 0.5, 0
    0         0    0.5
    0         0.85      0];


    shift  = 0;
    eg_rawtrace = figure;
    for i = ICuse
        plot(x,tc_avg(:,i)+shift,'Color',colord(mod(i-1,size(colord,1))+1,:),'linewidth', 1); hold on;
        %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
        %         hline(shift);
        %     end
        shift = shift+max(tc_avg(:,i))+200;  %10000 for 2P
    end
    set(gca,'YTick',5500:5500:shift);
    set(gca,'YTicklabel',(ICuse));
    
    xlabel('time(s)');ylabel('cell #');
    xlim([min(x) max(x)]);
    a = get(gca,'XTickLabel');
    set(gca,'XTickLabel',a,'FontSize',7);
    eg_rawtrace.Units = 'centimeters';
    eg_rawtrace.Position = [1 3 7 8];
    fig_name = ['eg_raw_traces_190603_img1025_20300-21100'];
    path = 'Z:\Analysis\figures\figure2_2Prun\';
    orient(eg_rawtrace,'landscape')
    print(eg_rawtrace,[path,fig_name],'-r600','-depsc');
end