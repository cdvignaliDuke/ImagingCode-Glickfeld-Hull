clear all
file_info
out_base = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis\');
for id = [4 8 9 10]%1:size(mouseID,2)
    for rID  = 1:2
        dest_sub  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[date{id}, '_', runID{rID}, '_', mouseID{id}],'\');
        if exist(dest_sub)
            load([dest_sub '_spont_events.mat']);
            load([dest_sub '_evoked_events.mat']);
            
            nCells = length(events);
            nDraws = 15;
            sampleN = randperm(nCells, nDraws);
            figure;
            for ic = 1:nDraws
                subplot(3,5,ic)
                
                fchunk = events(sampleN(ic)).df_chunk;
                for ip = 1:size(fchunk,1)
                    plot(fchunk(ip,:))
                    hold on
                end
            end
            supertitle(['Mouse ', mouseID{id},'Cell ',num2str(ic)]);
        end
    end
end

