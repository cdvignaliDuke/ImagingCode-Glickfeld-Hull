clear
file_info_CRP;
doPlotexample = 0; % if flaggged, will plot single cell Time course and segmentation
doCombMask = 0; % if flagged, it'll try find the folder with the combined session results
for sub =  [3, 6, 10, 25, 27, 32, 21, 8, 15, 28, 30, 34]
    %[2,5, 9, 22, 23, 31, 3, 6, 10, 25, 27, 32] for day1 and dayN
    %1:length(dates) 
    %[2, 5, 9, 22, 23, 31] day1, [3, 6, 10, 25, 27, 32] dayN
    for rID = 1:2
        file_info_CRP;
        
        dateID = dates;
        
        behav_dir = 'Z:\home\andrew\Behavior\Data\';
        
        data_dir = fullfile('Z:\home\jake\Analysis\Cue_reward_pairing_analysis\2P\');
        
        if doCombMask == 1
            subfold = fullfile(data_dir,[days_pair_folder{sub} '\']);
            dates = days_post;
            mouse = days_1_mouse{sub}
        else
            subfold = fullfile(data_dir,[dateID{sub} '_' runID{rID} '_' mouseID{sub} '\']);
            mouse = mouseID{sub}
        end
        
        if exist(subfold)
            load([subfold, 'ROI_TCs.mat']);
            load([subfold, '_cue_movies.mat'])
            
            
            date = dates{sub};
            
            dataName   = dir([behav_dir, '*', '9', mouse(end-1:end), '-', date(1:6), '*']);
            load([behav_dir, dataName(end).name]);
            tc_dir  = subfold;
            dest = tc_dir;
            dest_sub = dest; 
            
            
            if doPlotexample == 0
%                 getTC_events;  %umcomment to run, extract time course for each cell for event of interest
                CuePair_2P_TC_quantification; % quantify cell types and event types
%                 getTC_Spike_CRP; %umcomment to run, extract spikes
%                 spike_quantification_CRP; %umcomment to run, quantify spikes for rates and PSTH
                
            else
                %%% plot example segmentation
                figure;
                %                 subplot(2,1,1);
                %                 ax1=imshow(mat2gray(mean(data(:,:,1:5000),3)));
                
                %                 ax1 = subplot(2,1,1);
                imshow(mat2gray(mean(img_reg(:,:,10000:11000),3))); %10000:11000
                truesize
                %                 colormap(ax1,'gray');
                set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
                title('mean registered img');
                
                %                 subplot(2,1,2);
                %                 imagesc((mean(img_reg(:,:,5000:15000),3)));
                %                 colormap(ax1,'gray');
                figure;
                cmax = 3;
                cmin = 0;
                %                 ax2=subplot(2,1,2);
                mask_final_temp = zeros(size(mask_final));
                mask_final_temp(find(mask_final>0)) = 1;
                
                mask_ff = reshape(mask_final, [sz(1) sz(2)]);
                mask_ff(mask_ff == 11) = -1;
                mask_final_temp(find(mask_ff==-1)) = 2.5;
                imagesq(reshape(mask_final_temp, [sz(1) sz(2)]));truesize
                clim([cmin cmax]);
                %                 set(ax2, 'Clim', [cmin cmax]);
                title('segmentation');
                set(gca, 'XTickLabel', '', 'YTickLabel', '', 'Xtick', 0,'Ytick',0)
                
                
                data_tc_spont = tc_avg;
                nIC = size(tc_avg,2);
                %find all events to get rate
                events = [];
                events_diff_x = zeros(size(data_tc_spont,1)-1,nIC);
                % events_diff_x2 = zeros(size(data_tc_spont,1)-2,nIC);
                events_ind = {};
                events_rate = zeros(1,nIC);
                
                for ic = 1:nIC
                    events_diff_x(:,ic) = diff(data_tc_spont(:,ic),[],1);
                    %     events_diff_x1(:,ic) = diff(data_tc_spont(:,ic),[],1);
                    %     events_diff_x2(:,ic) = diff(events_diff_x1(:,ic),[],1);
                end
                
                
                thresh = 2.1;
                
                normalization = 1;
                deconvtau = 0;
                dt = 1;
                
                for ic = 1:nIC
                    %     events_ind{ic} = find(events_diff_x(:,ic)>thresh(ic));
                    [~, events_ind{ic}, ~] = CellsortFindspikes(events_diff_x(:,ic), thresh, dt, deconvtau, normalization);
                    events_ind{ic}(find(diff(events_ind{ic})==1)+1)=[];
                    events_rate(1,ic) = length(events_ind{ic})./((double(ifi)./1000)*size(data_tc_spont,1));
                end
                nDraws = size(tc_avg,2);
                %                 sampleN = randperm(nIC, nDraws);
                sampleN = 1:nDraws;
                
                for ic = 11:15
                    %     subplot(2,2,ic)
                    tt = (1:size(data_tc_spont(:,sampleN(ic)),1)) * 33;
                    figure;
                    plot(tt, data_tc_spont(:,sampleN(ic)) + 6000, 'k');
                    %     plot( data_tc_spont(:,sampleN(ic)) + 2, 'k');
                    hold on;
                    plot(tt(1:end-1), events_diff_x(:,sampleN(ic)) - 3, 'b');
                    %     plot( events_diff_x1(:,sampleN(ic)), 'b');
                    scale = max(events_diff_x(:,sampleN(ic)));
                    
                    plot( events_ind{sampleN(ic)}*33, scale*0.7, 'k.');
                    %     plot( events_ind{sampleN(ic)}, scale*0.6, 'k.');
                    
                    
                end
                % thresh = 2.5;
                title('Raw Trace--Black  1st Derivative--Blue');
            end
        end
        close all
        clearvars -except sub doPlotexample doCombMask
        
    end
end