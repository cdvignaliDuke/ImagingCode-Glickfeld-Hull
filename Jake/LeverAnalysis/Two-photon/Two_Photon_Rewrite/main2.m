clear
file_info;
usFacs = 100;
%14
for sub = 1:20
    for rID = 1:2
        file_info;
        mouse = mouseID{sub}
        dateID = date;
        date = dateID{sub}
        data_dir = fullfile('Z:\home\jake\Analysis\2P Analysis\Lever');
        subfold = fullfile(data_dir,[dateID{sub} '_' mouseID{sub} '_run' runID{rID} '\']);
        % foldName = subfold.name;
        % data_dir = fullfile(data_dir,foldName,'\');
        dataName   = dir(fullfile(subfold,'*frame_times.mat'));
        tc_dir  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[dateID{sub}, '_', runID{rID}, '_', mouseID{sub}],'\');
        dest = tc_dir;
        dest_sub = dest;
       % dest_sub = Z:\home\jake\Analysis\2P Analysis\Ziye_2P_figure
        
        if exist(dest_sub)
            load([tc_dir, 'ROI_TCs.mat']);
%             load([tc_dir, '_npSubTCs.mat']);
            
            if exist(subfold)
                load([subfold, dataName.name]);
            elseif strcmp(mouse, 'img30') && strcmp(date, '151018')
                dataName = 'data-i930-151018-1846.mat';
                subfold = 'Z:\home\andrew\Behavior\Data\';
                %                 dataName   = dir(fullfile(subfold, ['*' behavID{sub} '*' date]));
                load([subfold, dataName]);
            elseif strcmp(mouse, 'img32') && strcmp(date, '151216')
                dataName = 'data-i932-151216-1746.mat';
                subfold = 'Z:\home\andrew\Behavior\Data\';
                %                 dataName   = dir(fullfile(subfold, ['*' behavID{sub} '*' date]));
                load([subfold, dataName]);
            elseif strcmp(mouse, 'img36') && strcmp(date, '160203')
                dataName = 'data-i936-160203-1535.mat';
                subfold = 'Z:\home\andrew\Behavior\Data\';
                %                 dataName   = dir(fullfile(subfold, ['*' behavID{sub} '*' date]));
                load([subfold, dataName]);
            elseif strcmp(mouse, 'img36') && strcmp(date, '160205')
                dataName = 'data-i936-160205-1506.mat';
                subfold = 'Z:\home\andrew\Behavior\Data\';
                %                 dataName   = dir(fullfile(subfold, ['*' behavID{sub} '*' date]));
                load([subfold, dataName]);
             elseif strcmp(mouse, 'img59') && strcmp(date, '161102')
                dataName = 'data-i959-161102-1604.mat';
                subfold = 'Z:\home\andrew\Behavior\Data\';
                %                 dataName   = dir(fullfile(subfold, ['*' behavID{sub} '*' date]));
                load([subfold, dataName]);
             elseif strcmp(mouse, 'img59') && strcmp(date, '161103')
                dataName = 'data-i959-161103-1719.mat';
                subfold = 'Z:\home\andrew\Behavior\Data\';
                %                 dataName   = dir(fullfile(subfold, ['*' behavID{sub} '*' date]));
                load([subfold, dataName]);
                
            else
                error('subject %s_%s and subNum is %d input not found\n', date, mouse, sub);
            end
%             
%                         getTC_events;
%                         TC_quantification;
%                          close all

%             getTC_Spike
%             spike_outcomes
             spike_quantification
            close all
            clearvars -except mouseID sub rID
        end
    end
end