function [peak, peak_base, resp_avg, resp_sem] = getPeak(movie, window, frame, ifi, default_p)
peak = zeros(size(movie,1), size(movie,2));
peak_base = zeros(size(movie,1), size(movie,2));
peak_amp = zeros(size(movie,1), size(movie,2));
for i = 1:size(movie,1)
    for j = 1:size(movie,2) %each cell
        tc = squeeze(movie(i, j, window));
        tc_tot = squeeze(movie(i,j,:));
        [maxx, indx] = max(abs(tc));
        indx = window(indx);
        if indx == min(window) || indx == max(window)
            indx = median(window);
        end
%         if indx == length(tc) || indx == 1
%             indx = median(window);
%         end
        peak(i,j) = mean(tc_tot( (indx - floor(50./double(ifi))) : (indx + floor(50./double(ifi)))) );
%         base_window = (frame + ((indx - floor(50./double(ifi)) - 1 -1)) - 6):(frame + ((indx - floor(50./double(ifi)) - 1 -1)) - 3);
%         base_end = indx - round(300./double(ifi));
%         base_end = indx - round(300./double(ifi));
         base_window = 1 : (frame - round(300./double(ifi)));
        peak_base(i,j) = squeeze(mean(movie(i,j,base_window),3));
        peak_amp(i,j) = peak(i,j) - peak_base(i,j);
    end
end
resp_avg = mean(peak_amp,1);
resp_sem = std(peak_amp,[],1)./sqrt(size(peak_amp,1));

end