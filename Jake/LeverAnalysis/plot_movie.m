function plot_movie(movie, sz, rm_baseline,ts )
% --- plot sequence of move frames
%  movie could be a vercor or a matrix of frames 

all_frames = [];

dim3 = size(movie,3);


if(dim3 ==1) % single trial
    tmp_movie =[];
    tmp_movie(1,:,:) = movie;
    movie = tmp_movie;
end
d1 = size(movie,1);
d2 = size(movie,3);
counter = 0;
for j=1:d1
    for i=1:d2
        counter = counter+1;
        c_frame= reshape(movie(j,:,i),sz);
        
        if(exist('rm_baseline', 'var') && rm_baseline)  % remove baseline separetly for each trial
            use_inx = find(ts<-100); % use activity from 100 ms before event 
            if(isempty(use_inx))
                use_inx =1; 
            end
            base_avg = mean(movie(j,:,use_inx),3);
            c_frame= c_frame - reshape( base_avg ,sz);
        end
        c_frame(:,end) = NaN; % to make a separation line between frames
        all_frames(:,:,1,counter) = c_frame;
    end
end

% make all plots in the same color scale
min_val = min(all_frames(:));
max_val = max(all_frames(:)<inf);  %added the condition that it must be less that inf  JH
% ---- plot 
montage(all_frames, 'Size', [d1 d2]);
colormap jet;
caxis([min_val max_val]);
colorbar;
axis on;

% --- set xticks
xl = xlim;
yl = ylim;
x_tick = linspace(xl(1)+sz(2)/2, xl(2)-sz(2)/2, length(ts));
set(gca, 'xtick',x_tick, 'xticklabel',[]);
lb = strread(num2str(ts),'%s');
for i=1:length(lb)
    text(x_tick(i), yl(end)+15, lb{i}, 'HorizontalAlignment', 'center');
end

% ----- set yticks to the trial number if more than one trial 
if(size(movie,1)>1);
    y_tick = linspace(yl(1)+sz(1)/2, yl(2)-sz(1)/2, size(movie,1));
    set(gca, 'ytick',y_tick, 'yticklabel',1:size(movie,1));
    ylabel('#Trial');
else
    set(gca, 'ytick',[], 'yticklabel',[]);
end
