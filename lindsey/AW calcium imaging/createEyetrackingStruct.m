function mouse = createEyetrackingStruct(doPlot);
    close all
    AWEyeDatasets_AW
    av = behavParamsAV;
    rc = behavConstsAV;
    min_hold = 2000;
    pre_event_time = 1000;
    post_release_time = 1500;
    post_target_time = 4000;
    trans_win_time = [150 650];
    mice = unique({expt.SubNum});
    nMice = length(mice);
    str = unique({expt.SubNum});
    values = cell2mat(cellfun(@str2num,str,'UniformOutput',false));
    mouse_str = ['i' strjoin(str,'_i')];
    mouse_ind = find(intersect(cell2mat({av.mouse}),values));
    mouse_str = [];
    s = zeros(1,nMice);
    for imouse = 1:nMice
        mouse_str = [mouse_str 'i' num2str(av(mouse_ind(imouse)).mouse) '_'];  
    end
    for iexp = 1:size(expt,2)
        disp(num2str(iexp))
        SubNum = expt(iexp).SubNum;
        date_name = expt(iexp).date;
        runs = expt(iexp).runs;
        time_mat = expt(iexp).time_mat;
        mouse_name = expt(iexp).mouse;
        eyeFolder = expt(iexp).folder;
        frame_rate = expt(iexp).frame_rate;
        nrun = size(runs,1);
        VisonB2 = expt(iexp).VisonB2;

        if VisonB2 == 1
% %             str1 = sprintf('%f,', av.mouse);
%             values = textscan(str1, '%f', 'delimiter', ',', 'EmptyValue', NaN);
            imouse = find(values == str2num(SubNum));
            s(:,imouse) = s(:,imouse)+1;
            prepush_frames = ceil(pre_event_time*(frame_rate/1000));
%             trans_win_frames = ceil(pre_event_time+trans_win_time*(frame_rate/1000));
            postpush_frames = ceil(min_hold*(frame_rate/1000));
            prerelease_frames = ceil(pre_event_time*(frame_rate/1000));
            postrelease_frames = ceil(post_release_time*(frame_rate/1000));
            pretarget_frames = ceil(pre_event_time*(frame_rate/1000));
            posttarget_frames = ceil(post_target_time*(frame_rate/1000));
        %% load and combine mworks files
            for irun = 1:nrun
                time = time_mat(irun,:);
                fn_mworks = fullfile(rc.pathStr, ['data-i' SubNum '-' date_name '-' time '.mat']);
                if irun == 1
                    input = mwLoadData(fn_mworks, [], []);
                else
                    try
                        input = [input mwLoadData(fn_mworks, [], [])];
                    catch
                        input2 = mwLoadData(fn_mworks, [], []);
                        inpNames1 = fieldnames(input);
                        inpNames2 = fieldnames(input2);
                        inpLong = gt(length(inpNames1),length(inpNames2));
                        if inpLong == 1
                            inpPlusInd = ismember(inpNames1,inpNames2);
                            inpPlus = inpNames1(~inpPlusInd);
                            for i = 1:length(inpPlus)
                                input2.(genvarname(inpPlus(i))) = cell(1,80);
                            end
                        else
                            inpPlusInd = ismember(inpNames2,inpNames1);
                            inpPlus = inpNames2(~inpPlusInd);
                            for i = 1:length(inpPlus)
                                input.(char(genvarname(inpPlus(i)))) = cell(1,80);
                            end
                        end
                        input = [input input2];
                    end
                end
            end
            input = concatenateDataBlocks(input);

            runstr = runs(1,:);
            if nrun>1
                for irun = 2:nrun
                    runstr = [runstr '-' runs(irun,:)];
                end
            end
            if rc.name == 'ashley'
                fnout = fullfile(rc.eyeOutputDir,mouse_name,'eye tracking',date_name,[mouse_name '-' date_name '-' runstr]);
                fnin = fullfile(rc.eyeInputDir,mouse_name,'eye tracking',date_name,[mouse_name '-' date_name '-' runstr]);
            else
                fnout = fullfile(rc.eyeOutputDir,'eye tracking',date_name,[mouse_name '-' date_name], [mouse_name '-' date_name '-' runstr]);
                fnin = fnout;
            end
            load([fnin '_pupil.mat']);


            %% plot eye traces align to press
            set(0,'defaultfigurepaperorientation','portrait');
            set(0,'defaultfigurepapersize',[8.5 11]);
            set(0,'defaultfigurepaperposition',[.25 .25 [8.5 11]-0.5]);

            holdIx = find(cell2mat(input.holdTimesMs)>min_hold);
            rad_mat_down_base = bsxfun(@rdivide, rad_mat_down, mean(rad_mat_down(1:15,:),1));
            centroid_mat_down_base = bsxfun(@minus, centroid_mat_down, mean(centroid_mat_down(1:15,:,:),1));
            centroid_mat_down_base(:,2,:) = -1*centroid_mat_down_base(:,2,:);
            %plot change in eye area align to press
            if doPlot
                figure;
                tt = (1-prepush_frames:postpush_frames)*(1000/frame_rate);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_down(:,holdIx),2), nanstd(rad_mat_down(:,holdIx),[],2)./sqrt(length(holdIx)));
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_down(1,holdIx),2)).*[0.8 1.1]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,holdIx),2), nanstd(rad_mat_down_base(:,holdIx),[],2)./sqrt(length(holdIx)));
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.1])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,holdIx),3)), squeeze(nanstd(centroid_mat_down(:,1,holdIx),[],3))./sqrt(length(holdIx)));
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_down(1,1,holdIx),3))).*[0.9 1.1]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,holdIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,holdIx),[],3))./sqrt(length(holdIx)));
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,holdIx),3)), squeeze(nanstd(centroid_mat_down(:,2,holdIx),[],3))./sqrt(length(holdIx)));
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_down(1,2,holdIx),3))).*[0.9 1.1]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,holdIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,holdIx),[],3))./sqrt(length(holdIx)));
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle('Align to lever down')
                print([fnout '_avg_pressalign.pdf'], '-dpdf');
            end
            %plot change in Pupil radius by block type- success only
            b1Ix = intersect(holdIx, intersect(find(strcmp(input.trialOutcomeCell,'success')), find(cell2mat(input.tBlock2TrialNumber)==0)));
            b2Ix = intersect(holdIx, intersect(find(strcmp(input.trialOutcomeCell,'success')), find(cell2mat(input.tBlock2TrialNumber))));

            if doPlot
                figure;
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_down(:,b1Ix),2), nanstd(rad_mat_down(:,b1Ix),[],2)./sqrt(length(b1Ix)), '-g');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_down(:,b2Ix),2), nanstd(rad_mat_down(:,b2Ix),[],2)./sqrt(length(b2Ix)), '-k');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_down(1,[b1Ix b2Ix]),2)).*[0.8 1.1]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,b1Ix),2), nanstd(rad_mat_down_base(:,b1Ix),[],2)./sqrt(length(b1Ix)), '-g');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,b2Ix),2), nanstd(rad_mat_down_base(:,b2Ix),[],2)./sqrt(length(b2Ix)), '-k');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.9 1.1])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_down(:,1,b1Ix),[],3))./sqrt((length(b1Ix))), '-g');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_down(:,1,b2Ix),[],3))./sqrt((length(b2Ix))), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_down(1,1,[b1Ix b2Ix]),3))).*[0.9 1.1]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_down_base(:,1,b1Ix),[],3))./sqrt((length(b1Ix))), '-g');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_down_base(:,1,b2Ix),[],3))./sqrt((length(b2Ix))), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_down(:,2,b1Ix),[],3))./sqrt((length(b1Ix))), '-g');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_down(:,2,b2Ix),[],3))./sqrt((length(b2Ix))), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_down(1,2,[b1Ix b2Ix]),3))).*[0.9 1.1]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_down_base(:,2,b1Ix),[],3))./sqrt((length(b1Ix))), '-g');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_down_base(:,2,b2Ix),[],3))./sqrt((length(b2Ix))), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle([mouse_name ' ' date_name '- Align to lever down- Black:auditory Green: visual'])
                print([fnout '_avg_pressalign_AV.pdf'], '-dpdf');
            end
            mouse(imouse).expt(s(:,imouse)).align(1).name = 'press';
            mouse(imouse).expt(s(:,imouse)).align(2).name = 'release';
            mouse(imouse).expt(s(:,imouse)).align(3).name = 'target';
            for i = 1:3
                mouse(imouse).expt(s(:,imouse)).align(i).av(1).name = 'visual';
                mouse(imouse).expt(s(:,imouse)).align(i).av(2).name = 'auditory';
                for ii = 1:2
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(1).name = 'all';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(2).name = 'hit';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(3).name = 'miss';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(4).name = 'FA';
                    mouse(imouse).expt(s(:,imouse)).align(i).av(ii).outcome(5).name = 'CR';
                end
            end
            rad_trans_win = prepush_frames+round(prepush_frames/2):prepush_frames*2;%trans_win_frames;%
            sust_win = size(rad_mat_down,1)-round(prepush_frames*.667):size(rad_mat_down,1);
            cen_trans_win = prepush_frames+round(prepush_frames/2):prepush_frames+round(prepush_frames/2);%trans_win_frames;
            pre_win = 1:prepush_frames;

            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_rad_pre = [nanmean(nanmean(rad_mat_down(pre_win,b1Ix),1),2) nanstd(nanmean(rad_mat_down(pre_win,b1Ix),1),[],2)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_rad_pre = [nanmean(nanmean(rad_mat_down(pre_win,b2Ix),1),2) nanstd(nanmean(rad_mat_down(pre_win,b2Ix),1),[],2)/sqrt(length(b2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_hor_pre = [nanmean(nanmean(centroid_mat_down(pre_win,1,b1Ix),1),3) nanstd(nanmean(centroid_mat_down(pre_win,1,b1Ix),1),[],3)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_hor_pre = [nanmean(nanmean(centroid_mat_down(pre_win,1,b2Ix),1),3) nanstd(nanmean(centroid_mat_down(pre_win,1,b2Ix),1),[],3)/sqrt(length(b2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_ver_pre = [nanmean(nanmean(centroid_mat_down(pre_win,2,b1Ix),1),3) nanstd(nanmean(centroid_mat_down(pre_win,2,b1Ix),1),[],3)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_ver_pre = [nanmean(nanmean(centroid_mat_down(pre_win,2,b2Ix),1),3) nanstd(nanmean(centroid_mat_down(pre_win,2,b2Ix),1),[],3)/sqrt(length(b2Ix))];

            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_rad_trans = [nanmean(nanmean(rad_mat_down(rad_trans_win,b1Ix),1),2) nanstd(nanmean(rad_mat_down(rad_trans_win,b1Ix),1),[],2)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_rad_trans = [nanmean(nanmean(rad_mat_down(rad_trans_win,b2Ix),1),2) nanstd(nanmean(rad_mat_down(rad_trans_win,b2Ix),1),[],2)/sqrt(length(b2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_hor_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,1,b1Ix),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,1,b1Ix),1),[],3)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_hor_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,1,b2Ix),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,1,b2Ix),1),[],3)/sqrt(length(b2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_ver_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,2,b1Ix),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,2,b1Ix),1),[],3)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_ver_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,2,b2Ix),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,2,b2Ix),1),[],3)/sqrt(length(b2Ix))];

            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_rad_sust = [nanmean(nanmean(rad_mat_down(sust_win,b1Ix),1),2) nanstd(nanmean(rad_mat_down(sust_win,b1Ix),1),[],2)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_rad_sust = [nanmean(nanmean(rad_mat_down(sust_win,b2Ix),1),2) nanstd(nanmean(rad_mat_down(sust_win),1),[],2)/sqrt(length(b2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_hor_sust = [nanmean(nanmean(centroid_mat_down(sust_win,1,b1Ix),1),3) nanstd(nanmean(centroid_mat_down(sust_win,1,b1Ix),1),[],3)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_hor_sust = [nanmean(nanmean(centroid_mat_down(sust_win,1,b2Ix),1),3) nanstd(nanmean(centroid_mat_down(sust_win,1,b2Ix),1),[],3)/sqrt(length(b2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(1).avg_ver_sust = [nanmean(nanmean(centroid_mat_down(sust_win,2,b1Ix),1),3) nanstd(nanmean(centroid_mat_down(sust_win,2,b1Ix),1),[],3)/sqrt(length(b1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(1).avg_ver_sust = [nanmean(nanmean(centroid_mat_down(sust_win,2,b2Ix),1),3) nanstd(nanmean(centroid_mat_down(sust_win,2,b2Ix),1),[],3)/sqrt(length(b2Ix))];

            %plot change in Pupil radius by outcome type
            successIx = find(strcmp(input.trialOutcomeCell,'success'));
            missedIx = find(strcmp(input.trialOutcomeCell,'ignore'));
            if doPlot
                figure;
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_down(:,successIx),2), nanstd(rad_mat_down(:,successIx),[],2)./sqrt(length(successIx)), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_down(:,missedIx),2), nanstd(rad_mat_down(:,missedIx),[],2)./sqrt(length(missedIx)), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_down(1,:),2)).*[0.8 1.1]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,successIx),2), nanstd(rad_mat_down_base(:,successIx),[],2)./sqrt(length(successIx)), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,missedIx),2), nanstd(rad_mat_down_base(:,missedIx),[],2)./sqrt(length(missedIx)), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.1])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down(:,1,successIx),[],3))./sqrt((length(successIx))), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,1,missedIx),[],3))./sqrt((length(missedIx))), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_down(1,1,:),3))).*[0.9 1.1]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,successIx),[],3))./sqrt((length(successIx))), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,missedIx),[],3))./sqrt((length(missedIx))), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down(:,2,successIx),[],3))./sqrt((length(successIx))), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,2,missedIx),[],3))./sqrt((length(missedIx))), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_down(1,2,:),3))).*[0.9 1.1]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,successIx),[],3))./sqrt((length(successIx))), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,missedIx),[],3))./sqrt((length(missedIx))), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle([mouse_name ' ' date_name '- Align to lever down- Black: success; Red: missed'])
                print([fnout '_avg_pressalign_SM.pdf'], '-dpdf');
            end

            %hit and miss for V trials only
            b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
            b2Ix = find(cell2mat(input.tBlock2TrialNumber));
            successIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'success')));
            missedIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'ignore')));
            failureIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'failure')));
            if doPlot
                figure;
                downTrS = sum(~isnan(rad_mat_down(1,successIx)),2);
                downTrM = sum(~isnan(rad_mat_down(1,missedIx)),2);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_down(:,successIx),2), nanstd(rad_mat_down(:,successIx),[],2)./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_down(:,missedIx),2), nanstd(rad_mat_down(:,missedIx),[],2)./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_down(1,:),2)).*[0.8 1.4]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,successIx),2), nanstd(rad_mat_down_base(:,successIx),[],2)./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,missedIx),2), nanstd(rad_mat_down_base(:,missedIx),[],2)./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.4])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down(:,1,successIx),[],3))./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,1,missedIx),[],3))./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_down(1,1,:),3))).*[0.8 1.4]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,successIx),[],3))./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,missedIx),[],3))./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down(:,2,successIx),[],3))./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,2,missedIx),[],3))./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_down(1,2,:),3))).*[0.8 1.4]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,successIx),[],3))./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,missedIx),[],3))./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle([mouse_name ' ' date_name '- Align to press- Visual only- Black: success; Red: missed'])
                print([fnout '_avg_pressalign_SM_Vonly.pdf'], '-dpdf');
            end

            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_rad_pre = [nanmean(nanmean(rad_mat_down(pre_win,successIx),1),2) nanstd(nanmean(rad_mat_down(pre_win,successIx),1),[],2)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_rad_pre = [nanmean(nanmean(rad_mat_down(pre_win,missedIx),1),2) nanstd(nanmean(rad_mat_down(pre_win,missedIx),1),[],2)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(4).avg_rad_pre = [nanmean(nanmean(rad_mat_down(pre_win,failureIx),1),2) nanstd(nanmean(rad_mat_down(pre_win,failureIx),1),[],2)/sqrt(length(failureIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_hor_pre = [nanmean(nanmean(centroid_mat_down(pre_win,1,successIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,1,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_hor_pre = [nanmean(nanmean(centroid_mat_down(pre_win,1,missedIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,1,missedIx),1),[],3)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(4).avg_hor_pre = [nanmean(nanmean(centroid_mat_down(pre_win,1,failureIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,1,failureIx),1),[],3)/sqrt(length(failureIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_ver_pre = [nanmean(nanmean(centroid_mat_down(pre_win,2,successIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,2,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_ver_pre = [nanmean(nanmean(centroid_mat_down(pre_win,2,missedIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,2,missedIx),1),[],3)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(4).avg_ver_pre = [nanmean(nanmean(centroid_mat_down(pre_win,2,failureIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,2,failureIx),1),[],3)/sqrt(length(failureIx))];

            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_rad_trans = [nanmean(nanmean(rad_mat_down(rad_trans_win,successIx),1),2) nanstd(nanmean(rad_mat_down(rad_trans_win,successIx),1),[],2)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_rad_trans = [nanmean(nanmean(rad_mat_down(rad_trans_win,missedIx),1),2) nanstd(nanmean(rad_mat_down(rad_trans_win,missedIx),1),[],2)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_hor_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,1,successIx),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,1,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_hor_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,1,missedIx),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,1,missedIx),1),[],3)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_ver_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,2,successIx),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,2,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_ver_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,2,missedIx),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,2,missedIx),1),[],3)/sqrt(length(missedIx))];

            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_rad_sust = [nanmean(nanmean(rad_mat_down(sust_win,successIx),1),2) nanstd(nanmean(rad_mat_down(sust_win,successIx),1),[],2)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_rad_sust = [nanmean(nanmean(rad_mat_down(sust_win,missedIx),1),2) nanstd(nanmean(rad_mat_down(sust_win),1),[],2)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_hor_sust = [nanmean(nanmean(centroid_mat_down(sust_win,1,successIx),1),3) nanstd(nanmean(centroid_mat_down(sust_win,1,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_hor_sust = [nanmean(nanmean(centroid_mat_down(sust_win,1,missedIx),1),3) nanstd(nanmean(centroid_mat_down(sust_win,1,missedIx),1),[],3)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(2).avg_ver_sust = [nanmean(nanmean(centroid_mat_down(sust_win,2,successIx),1),3) nanstd(nanmean(centroid_mat_down(sust_win,2,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(1).outcome(3).avg_ver_sust = [nanmean(nanmean(centroid_mat_down(sust_win,2,missedIx),1),3) nanstd(nanmean(centroid_mat_down(sust_win,2,missedIx),1),[],3)/sqrt(length(missedIx))];

            %hit and miss for A trials only
            successIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
            missedIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'ignore')));
            failureIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'failure')));
            if doPlot
                figure;
                downTrS = sum(~isnan(rad_mat_down(1,successIx)),2);
                downTrM = sum(~isnan(rad_mat_down(1,missedIx)),2);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_down(:,successIx),2), nanstd(rad_mat_down(:,successIx),[],2)./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_down(:,missedIx),2), nanstd(rad_mat_down(:,missedIx),[],2)./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_down(1,:),2)).*[0.8 1.4]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,successIx),2), nanstd(rad_mat_down_base(:,successIx),[],2)./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_down_base(:,missedIx),2), nanstd(rad_mat_down_base(:,missedIx),[],2)./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.4])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down(:,1,successIx),[],3))./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,1,missedIx),[],3))./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_down(1,1,:),3))).*[0.8 1.4]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,successIx),[],3))./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,1,missedIx),[],3))./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down(:,2,successIx),[],3))./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down(:,2,missedIx),[],3))./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_down(1,2,:),3))).*[0.8 1.4]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,successIx),[],3))./sqrt(downTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_down_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_down_base(:,2,missedIx),[],3))./sqrt(downTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle([mouse_name ' ' date_name '- Align to press- Auditory only- Black: success; Red: missed'])
                print([fnout '_avg_pressalign_SM_Aonly.pdf'], '-dpdf');
            end

            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_rad_pre = [nanmean(nanmean(rad_mat_down(pre_win,successIx),1),2) nanstd(nanmean(rad_mat_down(pre_win,successIx),1),[],2)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_rad_pre = [nanmean(nanmean(rad_mat_down(pre_win,missedIx),1),2) nanstd(nanmean(rad_mat_down(pre_win,missedIx),1),[],2)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(4).avg_rad_pre = [nanmean(nanmean(rad_mat_down(pre_win,failureIx),1),2) nanstd(nanmean(rad_mat_down(pre_win,failureIx),1),[],2)/sqrt(length(failureIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_hor_pre = [nanmean(nanmean(centroid_mat_down(pre_win,1,successIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,1,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_hor_pre = [nanmean(nanmean(centroid_mat_down(pre_win,1,missedIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,1,missedIx),1),[],3)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(4).avg_hor_pre = [nanmean(nanmean(centroid_mat_down(pre_win,1,failureIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,1,failureIx),1),[],3)/sqrt(length(failureIx))];            
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_ver_pre = [nanmean(nanmean(centroid_mat_down(pre_win,2,successIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,2,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_ver_pre = [nanmean(nanmean(centroid_mat_down(pre_win,2,missedIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,2,missedIx),1),[],3)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(4).avg_ver_pre = [nanmean(nanmean(centroid_mat_down(pre_win,2,failureIx),1),3) nanstd(nanmean(centroid_mat_down(pre_win,2,failureIx),1),[],3)/sqrt(length(failureIx))];

            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_rad_trans = [nanmean(nanmean(rad_mat_down(rad_trans_win,successIx),1),2) nanstd(nanmean(rad_mat_down(rad_trans_win,successIx),1),[],2)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_rad_trans = [nanmean(nanmean(rad_mat_down(rad_trans_win,missedIx),1),2) nanstd(nanmean(rad_mat_down(rad_trans_win,missedIx),1),[],2)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_hor_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,1,successIx),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,1,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_hor_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,1,missedIx),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,1,missedIx),1),[],3)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_ver_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,2,successIx),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,2,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_ver_trans = [nanmean(nanmean(centroid_mat_down(cen_trans_win,2,missedIx),1),3) nanstd(nanmean(centroid_mat_down(cen_trans_win,2,missedIx),1),[],3)/sqrt(length(missedIx))];

            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_rad_sust = [nanmean(nanmean(rad_mat_down(sust_win,successIx),1),2) nanstd(nanmean(rad_mat_down(sust_win,successIx),1),[],2)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_rad_sust = [nanmean(nanmean(rad_mat_down(sust_win,missedIx),1),2) nanstd(nanmean(rad_mat_down(sust_win),1),[],2)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_hor_sust = [nanmean(nanmean(centroid_mat_down(sust_win,1,successIx),1),3) nanstd(nanmean(centroid_mat_down(sust_win,1,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_hor_sust = [nanmean(nanmean(centroid_mat_down(sust_win,1,missedIx),1),3) nanstd(nanmean(centroid_mat_down(sust_win,1,missedIx),1),[],3)/sqrt(length(missedIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(2).avg_ver_sust = [nanmean(nanmean(centroid_mat_down(sust_win,2,successIx),1),3) nanstd(nanmean(centroid_mat_down(sust_win,2,successIx),1),[],3)/sqrt(length(successIx))];
            mouse(imouse).expt(s(:,imouse)).align(1).av(2).outcome(3).avg_ver_sust = [nanmean(nanmean(centroid_mat_down(sust_win,2,missedIx),1),3) nanstd(nanmean(centroid_mat_down(sust_win,2,missedIx),1),[],3)/sqrt(length(missedIx))];



            %% plot change in pupil radius locked to lever up
            rad_mat_up_base = bsxfun(@rdivide, rad_mat_up, mean(rad_mat_down(1:15,:),1));
            centroid_mat_up_base = bsxfun(@minus, centroid_mat_up, mean(centroid_mat_down(1:15,:,:),1));
            centroid_mat_up_base(:,2,:) = -1*centroid_mat_up_base(:,2,:);
            
            tt = (1-prerelease_frames:postrelease_frames).*(1000/frame_rate);

            FIx = find(strcmp(input.trialOutcomeCell, 'failure'));
            SIx = find(strcmp(input.trialOutcomeCell, 'success'));
            MIx = find(strcmp(input.trialOutcomeCell, 'ignore'));
            FIxlong = intersect(find(cell2mat(input.tCyclesOn)>3), FIx);
            SIxlong = intersect(find(cell2mat(input.tCyclesOn)>3), SIx);
            b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
            b2Ix = find(cell2mat(input.tBlock2TrialNumber)==1);
            Fb1Ix = intersect(b1Ix, FIx);
            Fb2Ix = intersect(b2Ix, FIx);
            Sb1Ix = intersect(b1Ix, SIx);
            Sb2Ix = intersect(b2Ix, SIx);

            if doPlot
                figure;
                subplot(2,2,1)
                shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Fb1Ix),2), nanstd(rad_mat_up_base(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'r');
                hold on
                shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Sb1Ix),2), nanstd(rad_mat_up_base(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k'); 
                hold on
                vline(0,'k')
                title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
                subplot(2,2,2)
                shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Fb2Ix),2), nanstd(rad_mat_up_base(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'r');
                hold on
                shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Sb2Ix),2), nanstd(rad_mat_up_base(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'k'); 
                hold on
                vline(0,'k')
                title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
                subplot(2,2,3)
                shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Fb1Ix),2), nanstd(rad_mat_up_base(:,Fb1Ix),[],2)/sqrt(length(Fb1Ix)), 'k');
                hold on
                shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Fb2Ix),2), nanstd(rad_mat_up_base(:,Fb2Ix),[],2)/sqrt(length(Fb2Ix)), 'g'); 
                hold on
                vline(0,'k')
                title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
                subplot(2,2,4)
                shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Sb1Ix),2), nanstd(rad_mat_up_base(:,Sb1Ix),[],2)/sqrt(length(Sb1Ix)), 'k');
                hold on
                shadedErrorBar(tt, nanmean(rad_mat_up_base(:,Sb2Ix),2), nanstd(rad_mat_up_base(:,Sb2Ix),[],2)/sqrt(length(Sb2Ix)), 'g'); 
                hold on
                vline(0,'k')
                alignYaxes
                title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
                suptitle([mouse_name ' ' date_name '- Align to release- Radius']) 
                print([fnout '_avg_releasealign_SM_AV_rad.pdf'], '-dpdf');
                
                figure;
                subplot(2,2,1)
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,1,Fb1Ix),3), nanstd(centroid_mat_up_base(:,1,Fb1Ix),[],3)/sqrt(length(Fb1Ix)), 'r');
                hold on
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,1,Sb1Ix),3), nanstd(centroid_mat_up_base(:,1,Sb1Ix),[],3)/sqrt(length(Sb1Ix)), 'k'); 
                hold on
                vline(0,'k')
                title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
                subplot(2,2,2)
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,1,Fb2Ix),3), nanstd(centroid_mat_up_base(:,1,Fb2Ix),[],3)/sqrt(length(Fb2Ix)), 'r');
                hold on
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,1,Sb2Ix),3), nanstd(centroid_mat_up_base(:,1,Sb2Ix),[],3)/sqrt(length(Sb2Ix)), 'k'); 
                hold on
                vline(0,'k')
                title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
                subplot(2,2,3)
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,1,Fb1Ix),3), nanstd(centroid_mat_up_base(:,1,Fb1Ix),[],3)/sqrt(length(Fb1Ix)), 'k');
                hold on
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,1,Fb2Ix),3), nanstd(centroid_mat_up_base(:,1,Fb2Ix),[],3)/sqrt(length(Fb2Ix)), 'g'); 
                hold on
                vline(0,'k')
                title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
                subplot(2,2,4)
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,1,Sb1Ix),3), nanstd(centroid_mat_up_base(:,1,Sb1Ix),[],3)/sqrt(length(Sb1Ix)), 'k');
                hold on
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,1,Sb2Ix),3), nanstd(centroid_mat_up_base(:,1,Sb2Ix),[],3)/sqrt(length(Sb2Ix)), 'g'); 
                hold on
                vline(0,'k')
                alignYaxes
                title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
                suptitle([mouse_name ' ' date_name '- Align to release- Horizontal Pos']) 
                print([fnout '_avg_releasealign_SM_AV_hor.pdf'], '-dpdf');
                
                figure;
                subplot(2,2,1)
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,2,Fb1Ix),3), nanstd(centroid_mat_up_base(:,2,Fb1Ix),[],3)/sqrt(length(Fb1Ix)), 'r');
                hold on
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,2,Sb1Ix),3), nanstd(centroid_mat_up_base(:,2,Sb1Ix),[],3)/sqrt(length(Sb1Ix)), 'k'); 
                hold on
                vline(0,'k')
                title(['Visual trials: ' num2str(length(Sb1Ix)) ' Successes; ' num2str(length(Fb1Ix)) ' Earlies']) 
                subplot(2,2,2)
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,2,Fb2Ix),3), nanstd(centroid_mat_up_base(:,2,Fb2Ix),[],3)/sqrt(length(Fb2Ix)), 'r');
                hold on
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,2,Sb2Ix),3), nanstd(centroid_mat_up_base(:,2,Sb2Ix),[],3)/sqrt(length(Sb2Ix)), 'k'); 
                hold on
                vline(0,'k')
                title(['Auditory trials: ' num2str(length(Sb2Ix)) ' Successes; ' num2str(length(Fb2Ix)) ' Earlies']) 
                subplot(2,2,3)
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,2,Fb1Ix),3), nanstd(centroid_mat_up_base(:,2,Fb1Ix),[],3)/sqrt(length(Fb1Ix)), 'k');
                hold on
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,2,Fb2Ix),3), nanstd(centroid_mat_up_base(:,2,Fb2Ix),[],3)/sqrt(length(Fb2Ix)), 'g'); 
                hold on
                vline(0,'k')
                title(['Early trials: ' num2str(length(Fb1Ix)) ' Visual; ' num2str(length(Fb2Ix)) ' Auditory']) 
                subplot(2,2,4)
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,2,Sb1Ix),3), nanstd(centroid_mat_up_base(:,2,Sb1Ix),[],3)/sqrt(length(Sb1Ix)), 'k');
                hold on
                shadedErrorBar(tt, nanmean(centroid_mat_up_base(:,2,Sb2Ix),3), nanstd(centroid_mat_up_base(:,2,Sb2Ix),[],3)/sqrt(length(Sb2Ix)), 'g'); 
                hold on
                vline(0,'k')
                alignYaxes
                title(['Success trials: ' num2str(length(Sb1Ix)) ' Visual; ' num2str(length(Sb2Ix)) ' Auditory']) 
                print([fnout '_avg_releasealign_SM_AV_ver.pdf'], '-dpdf');
                suptitle([mouse_name ' ' date_name '- Align to release- Vertical Pos']) 
            end
            
            trans_win = prerelease_frames:prerelease_frames+round(prerelease_frames/3);
            sust_win = size(rad_mat_up,1)-round(prerelease_frames*0.667):size(rad_mat_up,1);
            pre_win = 1:prerelease_frames;
            
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_rad_pre = [nanmean(nanmean(rad_mat_up(pre_win,[Sb1Ix Fb1Ix]),1),2) nanstd(nanmean(rad_mat_up(pre_win,[Sb1Ix Fb1Ix]),1),[],2)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_rad_pre = [nanmean(nanmean(rad_mat_up(pre_win,[Sb2Ix Fb2Ix]),1),2) nanstd(nanmean(rad_mat_up(pre_win,[Sb2Ix Fb2Ix]),1),[],2)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_rad_pre = [nanmean(nanmean(rad_mat_up(pre_win,Sb1Ix),1),2) nanstd(nanmean(rad_mat_up(pre_win,Sb1Ix),1),[],2)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_rad_pre = [nanmean(nanmean(rad_mat_up(pre_win,Sb2Ix),1),2) nanstd(nanmean(rad_mat_up(pre_win,Sb2Ix),1),[],2)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_rad_pre = [nanmean(nanmean(rad_mat_up(pre_win,Fb1Ix),1),2) nanstd(nanmean(rad_mat_up(pre_win,Fb1Ix),1),[],2)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_rad_pre = [nanmean(nanmean(rad_mat_up(pre_win,Fb2Ix),1),2) nanstd(nanmean(rad_mat_up(pre_win,Fb2Ix),1),[],2)/sqrt(length(Fb2Ix))];
           
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_hor_pre = [nanmean(nanmean(centroid_mat_up(pre_win,1,[Sb1Ix Fb1Ix]),1),3) nanstd(nanmean(centroid_mat_up(pre_win,1,[Sb1Ix Fb1Ix]),1),[],3)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_hor_pre = [nanmean(nanmean(centroid_mat_up(pre_win,1,[Sb2Ix Fb2Ix]),1),3) nanstd(nanmean(centroid_mat_up(pre_win,1,[Sb2Ix Fb2Ix]),1),[],3)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_ver_pre = [nanmean(nanmean(centroid_mat_up(pre_win,2,[Sb1Ix Fb1Ix]),1),3) nanstd(nanmean(centroid_mat_up(pre_win,2,[Sb1Ix Fb1Ix]),1),[],3)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_ver_pre = [nanmean(nanmean(centroid_mat_up(pre_win,2,[Sb2Ix Fb2Ix]),1),3) nanstd(nanmean(centroid_mat_up(pre_win,2,[Sb2Ix Fb2Ix]),1),[],3)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_hor_pre = [nanmean(nanmean(centroid_mat_up(pre_win,1,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_up(pre_win,1,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_hor_pre = [nanmean(nanmean(centroid_mat_up(pre_win,1,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_up(pre_win,1,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_ver_pre = [nanmean(nanmean(centroid_mat_up(pre_win,2,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_up(pre_win,2,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_ver_pre = [nanmean(nanmean(centroid_mat_up(pre_win,2,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_up(pre_win,2,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_hor_pre = [nanmean(nanmean(centroid_mat_up(pre_win,1,Fb1Ix),1),3) nanstd(nanmean(centroid_mat_up(pre_win,1,Fb1Ix),1),[],3)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_hor_pre = [nanmean(nanmean(centroid_mat_up(pre_win,1,Fb2Ix),1),3) nanstd(nanmean(centroid_mat_up(pre_win,1,Fb2Ix),1),[],3)/sqrt(length(Fb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_ver_pre = [nanmean(nanmean(centroid_mat_up(pre_win,2,Fb1Ix),1),3) nanstd(nanmean(centroid_mat_up(pre_win,2,Fb1Ix),1),[],3)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_ver_pre = [nanmean(nanmean(centroid_mat_up(pre_win,2,Fb2Ix),1),3) nanstd(nanmean(centroid_mat_up(pre_win,2,Fb2Ix),1),[],3)/sqrt(length(Fb2Ix))];

            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_rad_trans = [nanmean(nanmean(rad_mat_up(trans_win,[Sb1Ix Fb1Ix]),1),2) nanstd(nanmean(rad_mat_up(trans_win,[Sb1Ix Fb1Ix]),1),[],2)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_rad_trans = [nanmean(nanmean(rad_mat_up(trans_win,[Sb2Ix Fb2Ix]),1),2) nanstd(nanmean(rad_mat_up(trans_win,[Sb2Ix Fb2Ix]),1),[],2)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_rad_trans = [nanmean(nanmean(rad_mat_up(trans_win,Sb1Ix),1),2) nanstd(nanmean(rad_mat_up(trans_win,Sb1Ix),1),[],2)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_rad_trans = [nanmean(nanmean(rad_mat_up(trans_win,Sb2Ix),1),2) nanstd(nanmean(rad_mat_up(trans_win,Sb2Ix),1),[],2)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_rad_trans = [nanmean(nanmean(rad_mat_up(trans_win,Fb1Ix),1),2) nanstd(nanmean(rad_mat_up(trans_win,Fb1Ix),1),[],2)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_rad_trans = [nanmean(nanmean(rad_mat_up(trans_win,Fb2Ix),1),2) nanstd(nanmean(rad_mat_up(trans_win,Fb2Ix),1),[],2)/sqrt(length(Fb2Ix))];
           
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_hor_trans = [nanmean(nanmean(centroid_mat_up(trans_win,1,[Sb1Ix Fb1Ix]),1),3) nanstd(nanmean(centroid_mat_up(trans_win,1,[Sb1Ix Fb1Ix]),1),[],3)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_hor_trans = [nanmean(nanmean(centroid_mat_up(trans_win,1,[Sb2Ix Fb2Ix]),1),3) nanstd(nanmean(centroid_mat_up(trans_win,1,[Sb2Ix Fb2Ix]),1),[],3)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_ver_trans = [nanmean(nanmean(centroid_mat_up(trans_win,2,[Sb1Ix Fb1Ix]),1),3) nanstd(nanmean(centroid_mat_up(trans_win,2,[Sb1Ix Fb1Ix]),1),[],3)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_ver_trans = [nanmean(nanmean(centroid_mat_up(trans_win,2,[Sb2Ix Fb2Ix]),1),3) nanstd(nanmean(centroid_mat_up(trans_win,2,[Sb2Ix Fb2Ix]),1),[],3)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_hor_trans = [nanmean(nanmean(centroid_mat_up(trans_win,1,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_up(trans_win,1,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_hor_trans = [nanmean(nanmean(centroid_mat_up(trans_win,1,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_up(trans_win,1,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_ver_trans = [nanmean(nanmean(centroid_mat_up(trans_win,2,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_up(trans_win,2,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_ver_trans = [nanmean(nanmean(centroid_mat_up(trans_win,2,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_up(trans_win,2,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_hor_trans = [nanmean(nanmean(centroid_mat_up(trans_win,1,Fb1Ix),1),3) nanstd(nanmean(centroid_mat_up(trans_win,1,Fb1Ix),1),[],3)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_hor_trans = [nanmean(nanmean(centroid_mat_up(trans_win,1,Fb2Ix),1),3) nanstd(nanmean(centroid_mat_up(trans_win,1,Fb2Ix),1),[],3)/sqrt(length(Fb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_ver_trans = [nanmean(nanmean(centroid_mat_up(trans_win,2,Fb1Ix),1),3) nanstd(nanmean(centroid_mat_up(trans_win,2,Fb1Ix),1),[],3)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_ver_trans = [nanmean(nanmean(centroid_mat_up(trans_win,2,Fb2Ix),1),3) nanstd(nanmean(centroid_mat_up(trans_win,2,Fb2Ix),1),[],3)/sqrt(length(Fb2Ix))];
    
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_rad_sust = [nanmean(nanmean(rad_mat_up(sust_win,[Sb1Ix Fb1Ix]),1),2) nanstd(nanmean(rad_mat_up(sust_win,[Sb1Ix Fb1Ix]),1),[],2)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_rad_sust = [nanmean(nanmean(rad_mat_up(sust_win,[Sb2Ix Fb2Ix]),1),2) nanstd(nanmean(rad_mat_up(sust_win,[Sb2Ix Fb2Ix]),1),[],2)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_rad_sust = [nanmean(nanmean(rad_mat_up(sust_win,Sb1Ix),1),2) nanstd(nanmean(rad_mat_up(sust_win,Sb1Ix),1),[],2)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_rad_sust = [nanmean(nanmean(rad_mat_up(sust_win,Sb2Ix),1),2) nanstd(nanmean(rad_mat_up(sust_win,Sb2Ix),1),[],2)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_rad_sust = [nanmean(nanmean(rad_mat_up(sust_win,Fb1Ix),1),2) nanstd(nanmean(rad_mat_up(sust_win,Fb1Ix),1),[],2)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_rad_sust = [nanmean(nanmean(rad_mat_up(sust_win,Fb2Ix),1),2) nanstd(nanmean(rad_mat_up(sust_win,Fb2Ix),1),[],2)/sqrt(length(Fb2Ix))];
           
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_hor_sust = [nanmean(nanmean(centroid_mat_up(sust_win,1,[Sb1Ix Fb1Ix]),1),3) nanstd(nanmean(centroid_mat_up(sust_win,1,[Sb1Ix Fb1Ix]),1),[],3)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_hor_sust = [nanmean(nanmean(centroid_mat_up(sust_win,1,[Sb2Ix Fb2Ix]),1),3) nanstd(nanmean(centroid_mat_up(sust_win,1,[Sb2Ix Fb2Ix]),1),[],3)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(1).avg_ver_sust = [nanmean(nanmean(centroid_mat_up(sust_win,2,[Sb1Ix Fb1Ix]),1),3) nanstd(nanmean(centroid_mat_up(sust_win,2,[Sb1Ix Fb1Ix]),1),[],3)/sqrt(length([Sb1Ix Fb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(1).avg_ver_sust = [nanmean(nanmean(centroid_mat_up(sust_win,2,[Sb2Ix Fb2Ix]),1),3) nanstd(nanmean(centroid_mat_up(sust_win,2,[Sb2Ix Fb2Ix]),1),[],3)/sqrt(length([Sb2Ix Fb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_hor_sust = [nanmean(nanmean(centroid_mat_up(sust_win,1,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_up(sust_win,1,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_hor_sust = [nanmean(nanmean(centroid_mat_up(sust_win,1,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_up(sust_win,1,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(2).avg_ver_sust = [nanmean(nanmean(centroid_mat_up(sust_win,2,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_up(sust_win,2,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(2).avg_ver_sust = [nanmean(nanmean(centroid_mat_up(sust_win,2,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_up(sust_win,2,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_hor_sust = [nanmean(nanmean(centroid_mat_up(sust_win,1,Fb1Ix),1),3) nanstd(nanmean(centroid_mat_up(sust_win,1,Fb1Ix),1),[],3)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_hor_sust = [nanmean(nanmean(centroid_mat_up(sust_win,1,Fb2Ix),1),3) nanstd(nanmean(centroid_mat_up(sust_win,1,Fb2Ix),1),[],3)/sqrt(length(Fb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(1).outcome(4).avg_ver_sust = [nanmean(nanmean(centroid_mat_up(sust_win,2,Fb1Ix),1),3) nanstd(nanmean(centroid_mat_up(sust_win,2,Fb1Ix),1),[],3)/sqrt(length(Fb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(2).av(2).outcome(4).avg_ver_sust = [nanmean(nanmean(centroid_mat_up(sust_win,2,Fb2Ix),1),3) nanstd(nanmean(centroid_mat_up(sust_win,2,Fb2Ix),1),[],3)/sqrt(length(Fb2Ix))];
    

        %% plot change in Pupil radius locked to target
            rad_mat_target_base = bsxfun(@rdivide, rad_mat_target, mean(rad_mat_down(1:15,:),1));
            centroid_mat_target_base = bsxfun(@minus, centroid_mat_target, mean(centroid_mat_down(1:15,:,:),1));
            centroid_mat_target_base(:,2,:) = -1*centroid_mat_target_base(:,2,:);
            
            if doPlot
                figure;
                tt = (1-pretarget_frames:posttarget_frames)*(1000/frame_rate);
                nonan_trials = sum(~isnan(rad_mat_target(1,:)),2);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_target,2), nanstd(rad_mat_target,[],2)./sqrt(nonan_trials));
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_target_base,2), nanstd(rad_mat_target_base,[],2)./sqrt(nonan_trials));
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.4])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,:),3)), squeeze(nanstd(centroid_mat_target(:,1,:),[],3))./sqrt(nonan_trials));
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,:),3)), squeeze(nanstd(centroid_mat_target_base(:,1,:),[],3))./sqrt(nonan_trials));
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,:),3)), squeeze(nanstd(centroid_mat_target(:,2,:),[],3))./sqrt(nonan_trials));
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                suptitle('Align to target')
                ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,:),3)), squeeze(nanstd(centroid_mat_target_base(:,2,:),[],3))./sqrt(nonan_trials));
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                print([fnout '_avg_targetalign.pdf'], '-dpdf');
            end

            %plot change in Pupil radius by block type
            b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
            b2Ix = find(cell2mat(input.tBlock2TrialNumber));

            if doPlot
                figure;
                targetTrB1 = sum(~isnan(rad_mat_target(1,b1Ix)),2);
                targetTrB2 = sum(~isnan(rad_mat_target(1,b2Ix)),2);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_target(:,b1Ix),2), nanstd(rad_mat_target(:,b1Ix),[],2)./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target(:,b2Ix),2), nanstd(rad_mat_target(:,b2Ix),[],2)./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b1Ix),2), nanstd(rad_mat_target_base(:,b1Ix),[],2)./sqrt(targetTrB1), '-g');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b2Ix),2), nanstd(rad_mat_target_base(:,b2Ix),[],2)./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.4])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b1Ix),[],3))./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b1Ix),[],3))./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b1Ix),[],3))./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b1Ix),[],3))./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle('Align to target- Black:auditory Green: visual')
                print([fnout '_avg_targetalign_AV.pdf'], '-dpdf');
            end

            %plot change in Pupil radius by block type- success only
            b1Ix = intersect(find(strcmp(input.trialOutcomeCell,'success')), find(cell2mat(input.tBlock2TrialNumber)==0));
            b2Ix = intersect(find(strcmp(input.trialOutcomeCell,'success')), find(cell2mat(input.tBlock2TrialNumber)));

            if doPlot
                figure;
                targetTrB1 = sum(~isnan(rad_mat_target(1,b1Ix)),2);
                targetTrB2 = sum(~isnan(rad_mat_target(1,b2Ix)),2);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_target(:,b1Ix),2), nanstd(rad_mat_target(:,b1Ix),[],2)./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target(:,b2Ix),2), nanstd(rad_mat_target(:,b2Ix),[],2)./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b1Ix),2), nanstd(rad_mat_target_base(:,b1Ix),[],2)./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,b2Ix),2), nanstd(rad_mat_target_base(:,b2Ix),[],2)./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.4])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b1Ix),[],3))./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b1Ix),[],3))./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,1,b2Ix),[],3))./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b1Ix),[],3))./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b1Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b1Ix),[],3))./sqrt(targetTrB1), 'g-');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,b2Ix),3)), squeeze(nanstd(centroid_mat_target_base(:,2,b2Ix),[],3))./sqrt(targetTrB2), '-k');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle('Align to target- Success only Black:auditory Green: visual')
                print([fnout '_avg_targetalign_AV_Sonly.pdf'], '-dpdf');
            end


            %plot change in Pupil radius by block type
            successIx = find(strcmp(input.trialOutcomeCell,'success'));
            missedIx = find(strcmp(input.trialOutcomeCell,'ignore'));

            if doPlot
                figure;
                targetTrS = sum(~isnan(rad_mat_target(1,successIx)),2);
                targetTrM = sum(~isnan(rad_mat_target(1,missedIx)),2);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_target(:,successIx),2), nanstd(rad_mat_target(:,successIx),[],2)./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target(:,missedIx),2), nanstd(rad_mat_target(:,missedIx),[],2)./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.4])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle('Align to target- Black: success; Red: missed')
                print([fnout '_avg_targetalign_SM.pdf'], '-dpdf');
            end

            %hit and miss for V trials only
            b1Ix = find(cell2mat(input.tBlock2TrialNumber)==0);
            b2Ix = find(cell2mat(input.tBlock2TrialNumber));
            successIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'success')));
            missedIx = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'ignore')));

            if doPlot
                figure;
                targetTrS = sum(~isnan(rad_mat_target(1,successIx)),2);
                targetTrM = sum(~isnan(rad_mat_target(1,missedIx)),2);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_target(:,successIx),2), nanstd(rad_mat_target(:,successIx),[],2)./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target(:,missedIx),2), nanstd(rad_mat_target(:,missedIx),[],2)./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.4])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle('Align to target- Visual only- Black: success; Red: missed')
                print([fnout '_avg_targetalign_SM_Vonly.pdf'], '-dpdf');
            end

            %hit and miss for A trials only
            successIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
            missedIx = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'ignore')));

            if doPlot
                figure;
                targetTrS = sum(~isnan(rad_mat_target(1,successIx)),2);
                targetTrM = sum(~isnan(rad_mat_target(1,missedIx)),2);
                subplot(3,2,1)
                shadedErrorBar(tt,nanmean(rad_mat_target(:,successIx),2), nanstd(rad_mat_target(:,successIx),[],2)./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target(:,missedIx),2), nanstd(rad_mat_target(:,missedIx),[],2)./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([(nanmean(rad_mat_target(1,:),2)).*[0.8 1.4]])
                subplot(3,2,2)
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,successIx),2), nanstd(rad_mat_target_base(:,successIx),[],2)./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,nanmean(rad_mat_target_base(:,missedIx),2), nanstd(rad_mat_target_base(:,missedIx),[],2)./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Pupil radius')
                ylim([0.8 1.4])
                subplot(3,2,3)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([squeeze((nanmean(centroid_mat_target(1,1,:),3))).*[0.8 1.4]])
                subplot(3,2,4)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,1,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,1,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Horizontal')
                ylim([-.1 .1])
                subplot(3,2,5)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([squeeze((nanmean(centroid_mat_target(1,2,:),3))).*[0.8 1.4]])
                subplot(3,2,6)
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,successIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,successIx),[],3))./sqrt(targetTrS), '-k');
                hold on
                shadedErrorBar(tt,squeeze(nanmean(centroid_mat_target_base(:,2,missedIx),3)), squeeze(nanstd(centroid_mat_target_base(:,2,missedIx),[],3))./sqrt(targetTrM), '-r');
                xlabel('Time (ms)')
                ylabel('Eye position- Vertical')
                ylim([-.1 .1])
                suptitle('Align to target- Auditory only- Black: success; Red: missed')
                print([fnout '_avg_targetalign_SM_Aonly.pdf'], '-dpdf');
            end
            
            trans_win = pretarget_frames+round(pretarget_frames/7):pretarget_frames+round(pretarget_frames/2);
            sust_win = 2*(pretarget_frames):3*(pretarget_frames);
            pre_win = 1:pretarget_frames;
            Sb1Ix = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'success')));
            Mb1Ix = intersect(b1Ix, find(strcmp(input.trialOutcomeCell,'ignore')));
            Sb2Ix = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'success')));
            Mb2Ix = intersect(b2Ix, find(strcmp(input.trialOutcomeCell,'ignore')));

            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_rad_pre = [nanmean(nanmean(rad_mat_target(pre_win,[Sb1Ix Mb1Ix]),1),2) nanstd(nanmean(rad_mat_target(pre_win,[Sb1Ix Mb1Ix]),1),[],2)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_rad_pre = [nanmean(nanmean(rad_mat_target(pre_win,[Sb2Ix Mb2Ix]),1),2) nanstd(nanmean(rad_mat_target(pre_win,[Sb2Ix Mb2Ix]),1),[],2)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_rad_pre = [nanmean(nanmean(rad_mat_target(pre_win,Sb1Ix),1),2) nanstd(nanmean(rad_mat_target(pre_win,Sb1Ix),1),[],2)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_rad_pre = [nanmean(nanmean(rad_mat_target(pre_win,Sb2Ix),1),2) nanstd(nanmean(rad_mat_target(pre_win,Sb2Ix),1),[],2)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_rad_pre = [nanmean(nanmean(rad_mat_target(pre_win,Mb1Ix),1),2) nanstd(nanmean(rad_mat_target(pre_win,Mb1Ix),1),[],2)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_rad_pre = [nanmean(nanmean(rad_mat_target(pre_win,Mb2Ix),1),2) nanstd(nanmean(rad_mat_target(pre_win,Mb2Ix),1),[],2)/sqrt(length(Mb2Ix))];
           
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_hor_pre = [nanmean(nanmean(centroid_mat_target(pre_win,1,[Sb1Ix Mb1Ix]),1),3) nanstd(nanmean(centroid_mat_target(pre_win,1,[Sb1Ix Mb1Ix]),1),[],3)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_hor_pre = [nanmean(nanmean(centroid_mat_target(pre_win,1,[Sb2Ix Mb2Ix]),1),3) nanstd(nanmean(centroid_mat_target(pre_win,1,[Sb2Ix Mb2Ix]),1),[],3)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_ver_pre = [nanmean(nanmean(centroid_mat_target(pre_win,2,[Sb1Ix Mb1Ix]),1),3) nanstd(nanmean(centroid_mat_target(pre_win,2,[Sb1Ix Mb1Ix]),1),[],3)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_ver_pre = [nanmean(nanmean(centroid_mat_target(pre_win,2,[Sb2Ix Mb2Ix]),1),3) nanstd(nanmean(centroid_mat_target(pre_win,2,[Sb2Ix Mb2Ix]),1),[],3)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_hor_pre = [nanmean(nanmean(centroid_mat_target(pre_win,1,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_target(pre_win,1,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_hor_pre = [nanmean(nanmean(centroid_mat_target(pre_win,1,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_target(pre_win,1,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_ver_pre = [nanmean(nanmean(centroid_mat_target(pre_win,2,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_target(pre_win,2,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_ver_pre = [nanmean(nanmean(centroid_mat_target(pre_win,2,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_target(pre_win,2,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_hor_pre = [nanmean(nanmean(centroid_mat_target(pre_win,1,Mb1Ix),1),3) nanstd(nanmean(centroid_mat_target(pre_win,1,Mb1Ix),1),[],3)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_hor_pre = [nanmean(nanmean(centroid_mat_target(pre_win,1,Mb2Ix),1),3) nanstd(nanmean(centroid_mat_target(pre_win,1,Mb2Ix),1),[],3)/sqrt(length(Mb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_ver_pre = [nanmean(nanmean(centroid_mat_target(pre_win,2,Mb1Ix),1),3) nanstd(nanmean(centroid_mat_target(pre_win,2,Mb1Ix),1),[],3)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_ver_pre = [nanmean(nanmean(centroid_mat_target(pre_win,2,Mb2Ix),1),3) nanstd(nanmean(centroid_mat_target(pre_win,2,Mb2Ix),1),[],3)/sqrt(length(Mb2Ix))];

            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_rad_trans = [nanmean(nanmean(rad_mat_target(trans_win,[Sb1Ix Mb1Ix]),1),2) nanstd(nanmean(rad_mat_target(trans_win,[Sb1Ix Mb1Ix]),1),[],2)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_rad_trans = [nanmean(nanmean(rad_mat_target(trans_win,[Sb2Ix Mb2Ix]),1),2) nanstd(nanmean(rad_mat_target(trans_win,[Sb2Ix Mb2Ix]),1),[],2)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_rad_trans = [nanmean(nanmean(rad_mat_target(trans_win,Sb1Ix),1),2) nanstd(nanmean(rad_mat_target(trans_win,Sb1Ix),1),[],2)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_rad_trans = [nanmean(nanmean(rad_mat_target(trans_win,Sb2Ix),1),2) nanstd(nanmean(rad_mat_target(trans_win,Sb2Ix),1),[],2)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_rad_trans = [nanmean(nanmean(rad_mat_target(trans_win,Mb1Ix),1),2) nanstd(nanmean(rad_mat_target(trans_win,Mb1Ix),1),[],2)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_rad_trans = [nanmean(nanmean(rad_mat_target(trans_win,Mb2Ix),1),2) nanstd(nanmean(rad_mat_target(trans_win,Mb2Ix),1),[],2)/sqrt(length(Mb2Ix))];
           
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_hor_trans = [nanmean(nanmean(centroid_mat_target(trans_win,1,[Sb1Ix Mb1Ix]),1),3) nanstd(nanmean(centroid_mat_target(trans_win,1,[Sb1Ix Mb1Ix]),1),[],3)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_hor_trans = [nanmean(nanmean(centroid_mat_target(trans_win,1,[Sb2Ix Mb2Ix]),1),3) nanstd(nanmean(centroid_mat_target(trans_win,1,[Sb2Ix Mb2Ix]),1),[],3)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_ver_trans = [nanmean(nanmean(centroid_mat_target(trans_win,2,[Sb1Ix Mb1Ix]),1),3) nanstd(nanmean(centroid_mat_target(trans_win,2,[Sb1Ix Mb1Ix]),1),[],3)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_ver_trans = [nanmean(nanmean(centroid_mat_target(trans_win,2,[Sb2Ix Mb2Ix]),1),3) nanstd(nanmean(centroid_mat_target(trans_win,2,[Sb2Ix Mb2Ix]),1),[],3)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_hor_trans = [nanmean(nanmean(centroid_mat_target(trans_win,1,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_target(trans_win,1,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_hor_trans = [nanmean(nanmean(centroid_mat_target(trans_win,1,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_target(trans_win,1,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_ver_trans = [nanmean(nanmean(centroid_mat_target(trans_win,2,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_target(trans_win,2,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_ver_trans = [nanmean(nanmean(centroid_mat_target(trans_win,2,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_target(trans_win,2,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_hor_trans = [nanmean(nanmean(centroid_mat_target(trans_win,1,Mb1Ix),1),3) nanstd(nanmean(centroid_mat_target(trans_win,1,Mb1Ix),1),[],3)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_hor_trans = [nanmean(nanmean(centroid_mat_target(trans_win,1,Mb2Ix),1),3) nanstd(nanmean(centroid_mat_target(trans_win,1,Mb2Ix),1),[],3)/sqrt(length(Mb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_ver_trans = [nanmean(nanmean(centroid_mat_target(trans_win,2,Mb1Ix),1),3) nanstd(nanmean(centroid_mat_target(trans_win,2,Mb1Ix),1),[],3)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_ver_trans = [nanmean(nanmean(centroid_mat_target(trans_win,2,Mb2Ix),1),3) nanstd(nanmean(centroid_mat_target(trans_win,2,Mb2Ix),1),[],3)/sqrt(length(Mb2Ix))];
    
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_rad_sust = [nanmean(nanmean(rad_mat_target(sust_win,[Sb1Ix Mb1Ix]),1),2) nanstd(nanmean(rad_mat_target(sust_win,[Sb1Ix Mb1Ix]),1),[],2)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_rad_sust = [nanmean(nanmean(rad_mat_target(sust_win,[Sb2Ix Mb2Ix]),1),2) nanstd(nanmean(rad_mat_target(sust_win,[Sb2Ix Mb2Ix]),1),[],2)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_rad_sust = [nanmean(nanmean(rad_mat_target(sust_win,Sb1Ix),1),2) nanstd(nanmean(rad_mat_target(sust_win,Sb1Ix),1),[],2)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_rad_sust = [nanmean(nanmean(rad_mat_target(sust_win,Sb2Ix),1),2) nanstd(nanmean(rad_mat_target(sust_win,Sb2Ix),1),[],2)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_rad_sust = [nanmean(nanmean(rad_mat_target(sust_win,Mb1Ix),1),2) nanstd(nanmean(rad_mat_target(sust_win,Mb1Ix),1),[],2)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_rad_sust = [nanmean(nanmean(rad_mat_target(sust_win,Mb2Ix),1),2) nanstd(nanmean(rad_mat_target(sust_win,Mb2Ix),1),[],2)/sqrt(length(Mb2Ix))];
           
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_hor_sust = [nanmean(nanmean(centroid_mat_target(sust_win,1,[Sb1Ix Mb1Ix]),1),3) nanstd(nanmean(centroid_mat_target(sust_win,1,[Sb1Ix Mb1Ix]),1),[],3)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_hor_sust = [nanmean(nanmean(centroid_mat_target(sust_win,1,[Sb2Ix Mb2Ix]),1),3) nanstd(nanmean(centroid_mat_target(sust_win,1,[Sb2Ix Mb2Ix]),1),[],3)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(1).avg_ver_sust = [nanmean(nanmean(centroid_mat_target(sust_win,2,[Sb1Ix Mb1Ix]),1),3) nanstd(nanmean(centroid_mat_target(sust_win,2,[Sb1Ix Mb1Ix]),1),[],3)/sqrt(length([Sb1Ix Mb1Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(1).avg_ver_sust = [nanmean(nanmean(centroid_mat_target(sust_win,2,[Sb2Ix Mb2Ix]),1),3) nanstd(nanmean(centroid_mat_target(sust_win,2,[Sb2Ix Mb2Ix]),1),[],3)/sqrt(length([Sb2Ix Mb2Ix]))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_hor_sust = [nanmean(nanmean(centroid_mat_target(sust_win,1,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_target(sust_win,1,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_hor_sust = [nanmean(nanmean(centroid_mat_target(sust_win,1,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_target(sust_win,1,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(2).avg_ver_sust = [nanmean(nanmean(centroid_mat_target(sust_win,2,Sb1Ix),1),3) nanstd(nanmean(centroid_mat_target(sust_win,2,Sb1Ix),1),[],3)/sqrt(length(Sb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(2).avg_ver_sust = [nanmean(nanmean(centroid_mat_target(sust_win,2,Sb2Ix),1),3) nanstd(nanmean(centroid_mat_target(sust_win,2,Sb2Ix),1),[],3)/sqrt(length(Sb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_hor_sust = [nanmean(nanmean(centroid_mat_target(sust_win,1,Mb1Ix),1),3) nanstd(nanmean(centroid_mat_target(sust_win,1,Mb1Ix),1),[],3)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_hor_sust = [nanmean(nanmean(centroid_mat_target(sust_win,1,Mb2Ix),1),3) nanstd(nanmean(centroid_mat_target(sust_win,1,Mb2Ix),1),[],3)/sqrt(length(Mb2Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(1).outcome(3).avg_ver_sust = [nanmean(nanmean(centroid_mat_target(sust_win,2,Mb1Ix),1),3) nanstd(nanmean(centroid_mat_target(sust_win,2,Mb1Ix),1),[],3)/sqrt(length(Mb1Ix))];
            mouse(imouse).expt(s(:,imouse)).align(3).av(2).outcome(3).avg_ver_sust = [nanmean(nanmean(centroid_mat_target(sust_win,2,Mb2Ix),1),3) nanstd(nanmean(centroid_mat_target(sust_win,2,Mb2Ix),1),[],3)/sqrt(length(Mb2Ix))]; 
            
       end
    end
%     mouse_str = [];
%     for imouse = 1:size(av,2)
%         mouse_str = [mouse_str 'i' num2str(av(imouse).mouse) '_'];  
%     end
   save(fullfile(rc.eyeOutputDir, [date '_' mouse_str 'EyeSummary.mat']), 'mouse');

end

