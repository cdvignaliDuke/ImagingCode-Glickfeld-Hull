load('cell_needTOchangeINTOmatrix.mat')
load('speed.mat')
vel=[]; % a line of velocity value 
flu=[]; % a line of fluorescence value 
for i=1:1:6
    a=dfOvF_spd_all{i};
    flu=[flu;a(:)];
    [la1,la2]=size(a);
    b=isp_mat(i,:);
    b(isnan(b))=[]; %delate NaN
    b2=repmat(b,la1,1);
    vel=[vel;b2(:)];    
end 
vel_plot_x=sort(unique(vel));       %reuslt x
flu_plot_y=zeros(length(vel_plot_x),1); %result y
errorbar_plot_y=zeros(length(vel_plot_x),1); %result y errorbar
for j=1:1:length(vel_plot_x)
    c=find(vel==vel_plot_x(j));
    flu_plot_y(j)=sum(flu(c))/length(c);
    errorbar_plot_y(j)=std(flu(c))/sqrt(length(c));
end 
plot(vel_plot_x,flu_plot_y)
