function p=tcPlotTcourses(expt, tcnorm, avg_img, labelimg, titleName, print_flag, ax)

% expt: structure
%   needs expt.trialdur, expt.ntrials, expt.stimdur
% tcnorm: normalized timecourses around 1
% avg_img: average stack image
% labelimg: cell labels
% titleName: (optional) whatever you want
% print_flag: (default=1) 1: printer, 2: postscript 3: both
%
% modified form SY's original code
% Kenichi Ohki  3/18/08

if nargin < 5
    titleName = 'NoName';
end
if nargin < 6
    print_flag = 1;
end
if nargin < 7    
    ax=[];
end

p=0;
debugflag=false;


sz=size(tcnorm);

tcavgnorm = squeeze(mean(reshape(tcnorm,expt.trialdur,expt.ntrials,sz(2)),2));

n_perpage=15;

npage=ceil(sz(2)/n_perpage);

n=1;

warning off all        

avg_img=double(avg_img);
ker=fspecial('gaussian',90,30);
avg_f=imFilter2(ker,avg_img);
avg_img_2=avg_img./avg_f;
avg_img_2=avg_img_2-min(min(avg_img_2));

msk=im2bw(labelimg,0.0001);
im3c=avg_img_2./max(max(avg_img_2));
im3c(:,:,2)=avg_img_2./max(max(avg_img_2));
im3c(:,:,3)=avg_img_2./max(max(avg_img_2));
im3c(:,:,1)=(avg_img_2./max(max(avg_img_2))).*(1-msk);

for i=1:npage
    h=figure;
    set(h,'PaperUnits','normalized');
    set(h,'PaperPosition',[0.05 0.05 0.9 0.9]);
    %main comment

    subplot('Position',[0.05 0.95 0.9 0.05]); 
    set(gca,'FontSize',6);
    axis off;

%    text(0.1,0.8,num2str(s.nline),'FontSize',6);
    text(0.1,0.6,[num2str(i) '/' num2str(npage)],'FontSize',6);
    text(0.1,0.4,titleName,'FontSize',6);
%    text(0.1,0.2,matoutfname,'FontSize',6);


    % image
    subplot('Position',[0.05 0.8 0.2 0.14]); 
    set(gca,'FontSize',6);
    %axis square;
    image(avg_img*64/max(max(avg_img)));
    colormap gray;

    % image2
    subplot('Position',[0.3 0.8 0.2 0.14]); 
    set(gca,'FontSize',6);
    %axis square;
    image(avg_img_2*64/max(max(avg_img_2)));
    image(im3c);

    %cells
    h=subplot('Position',[0.55 0.8 0.2 0.14]); 
    set(gca,'FontSize',6);
    %axis square

    map(1:sz(1),1:3)=0.5;
    map((i-1)*n_perpage+1:min(i*n_perpage,sz(2)),1:3)=1;
    lblim=label2rgb(labelimg,map,'k');
    image(lblim);

    hold on;

    for j=1:n_perpage
        if n>sz(2)
%                if print_flag==1
%                    print
%                end
%   this is not necessary now
%
            break;  % changed from return to break by KO 09/15/04
        end

        axes(h);
        [y,x] = find(labelimg==n);
        text(x(1),y(1),num2str(n),'Color','r','FontSize',6);
    
    
        %subcomment
        subplot('Position',[0.8 0.80-j*0.052 0.15 0.05]); 
        axis off;
        text(0.1,0.8,num2str(n));
    
        subplot('Position',[0.05 0.80-j*0.052 0.2 0.05]);
        set(gca,'FontSize',6);
        plot(tcnorm(:,n));

        set(gca,'xtick',0:expt.trialdur:sz(1));
        set(gca,'xgrid','on');
        if  ~isempty(ax)
            a=axis;
            a(3)=1-ax;
            a(4)=1+ax;
            axis(a);
        end
        %
        subplot('Position',[0.3 0.80-j*0.052 0.2 0.05]);
        set(gca,'FontSize',6);
        plot(tcavgnorm(:,n));
        set(gca,'xtick',0:expt.stimdur:length(tcavgnorm(:,n)));
        set(gca,'xgrid','on');

        if  ~isempty(ax)
            a=axis;
            a(3)=1-ax;
            a(4)=1+ax;
            axis(a);
        end
    
    
        %
        subplot('Position',[0.55 0.80-j*0.052 0.2 0.05]);
        set(gca,'FontSize',6);
        szlong=size(tcnorm);
        szshort=size(tcavgnorm);

        all_shorts=reshape(tcnorm(:,n),szshort(1),szlong(1)/szshort(1));
        plot(all_shorts);
        set(gca,'xtick',0:10:szshort(1));
        set(gca,'xgrid','on');

       if  ~isempty(ax)
            a=axis;
            a(3)=1-ax;
            a(4)=1+ax;
            axis(a);
        end
        n=n+1; 
    end 

    if print_flag==1
        print;
    end

    % output as a postscript file
    % added by KO 09/15/04

    if print_flag==2
        
        print ('-dpsc','-append',[titleName, '_tcourse_plot.ps']);
    end
    if print_flag==3
        print;
        print ('-dpsc','-append',[titleName, '_tcourse_plot.ps']);
    end

end
t=1;

warning on all   

if debugflag
    disp('Pause,  hit CR');
    pause
end
close all
end

  

