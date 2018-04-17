% get 10% rise time of averaged F over trials

function [riseIdx] = findRisetime(avg_movie, pre_frames)
if length(size(avg_movie)) < 3
    nCells = size(avg_movie,1);
    riseIdx = zeros(nCells,1);
    for i = 1:nCells % loop each cell
        
        %     [~, LT, ~] = risetime(avg_movie(i,pre_cue_frames:pre_cue_frames + 31));
        %     if ~isempty(LT)
        %         riseIdx(i) = round(LT(end));
        %     end
        
        movie_temp = avg_movie(i,pre_frames: pre_frames + 50);
        
        [maxv, riseIdx(i)] = max(movie_temp);
        
        k = riseIdx(i) - 1;
        
        while k >=1
            if movie_temp(k) > 0.25* maxv
                k = k - 1;
            else
                riseIdx(i) = k;
                k = 0;
            end
        end
    end
else
    nCells = size(avg_movie,2);
    nEvents = size(avg_movie,1);
    riseIdx = zeros(nEvents, nCells);
    for i = 1:nEvents
        for j = 1:nCells % loop each cell
            
            %     [~, LT, ~] = risetime(avg_movie(i,pre_cue_frames:pre_cue_frames + 31));
            %     if ~isempty(LT)
            %         riseIdx(i) = round(LT(end));
            %     end
            
            movie_temp = squeeze(avg_movie(i,j,:));
            
            [maxv, idx] = max(movie_temp);
            
            k = idx - 1;
            
            while k >=1
                if movie_temp(k) > 0.25* maxv
                    k = k - 1;
                else
                    riseIdx(i,j) = k - pre_frames;
                    k = 0;
                end
            end
        end
    end
end
% riseIdx(riseIdx < pre_cue_frames) = [];
end