function plot_movie(movie, sz, rm_baseline,ts )
% --- plot sequence of move frames
%  movie could be a vector or a matrix of frames 

%allocate memory, check size of movie and record dimensions
sz = sz(1:2);
all_frames = [];
dim3 = size(movie,3);
if(dim3 ==1) % single trial or the movie has already been averaged
    tmp_movie =[];
    tmp_movie(1,:,:) = movie;
    movie = tmp_movie;
end
dim1 = size(movie,1); %trial num
dim2 = size(movie,3); %frame num within trial

%concetenates all the frames into a rectangle with dividing lines (NaN)
counter = 0;
for j=1:dim1
    for ii=1:dim2
        counter = counter+1;
        c_frame= reshape(movie(j,:,ii),sz);
        if(exist('rm_baseline', 'var') && rm_baseline)  % remove baseline separetly for each trial
            use_inx = find(ts<-400); % use activity from 400 ms before event 
            if(isempty(use_inx))
                use_inx =1; 
            end
            base_avg = mean(movie(j,:,use_inx),3);
            c_frame= c_frame - reshape( base_avg ,sz);
        end
        c_frame(:,[end:end+1]) = NaN; % to make a separation line between frames
        all_frames(:,:,1,counter) = c_frame;
    end
end

% make all plots in the same color scale
min_val = min(all_frames(:));
max_val = max(all_frames(:)<inf);  %added the condition that it must be less that inf  JH
% ---- plot 
montage(all_frames, 'Size', [dim1 dim2]);
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
for ii=1:length(lb)
    text(x_tick(ii), yl(end)+15, lb{ii}, 'HorizontalAlignment', 'center');
end

% ----- set yticks to the trial number if more than one trial 
if(size(movie,1)>1);
    y_tick = linspace(yl(1)+sz(1)/2, yl(2)-sz(1)/2, size(movie,1));
    set(gca, 'ytick',y_tick, 'yticklabel',1:size(movie,1));
    ylabel('#Trial');
else
    set(gca, 'ytick',[], 'yticklabel',[]);
end
