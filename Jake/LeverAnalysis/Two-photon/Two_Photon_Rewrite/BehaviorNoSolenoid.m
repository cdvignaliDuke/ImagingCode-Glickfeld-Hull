file_info;
runID = {'000', '001', '001', '000','000'};
i = 2;
out_dir  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',[date{i}, '_', runID{i}, '_', mouseID{i}],'\');
data_dir = fullfile('\\crash\data\home\jake\Analysis\2P Analysis\Lever',[date{i} '_' mouseID{i} '_run' runID{i}]);

fileName   = dir(fullfile(data_dir,'*frame_times.mat'));

load([data_dir,'\',fileName.name]);

leverTimeC = input.leverTimesUs;

leverValueC = input.leverValues;

trialLen = length(leverValueC);
figure;
press_v = zeros(trialLen,1);
release_v = zeros(trialLen,1);

press_v_s = zeros(trialLen,1);
release_v_s = zeros(trialLen,1);

press_v_f = zeros(trialLen,1);
release_v_f = zeros(trialLen,1);
shift=0;
succ=0;
fail = 0;
% center timestamp to first release
for i =  1: trialLen
    if length(leverValueC{1,i}) > 2
        tempV = (leverValueC{1,i});
        
        tempT = (leverTimeC{1,i})/1000;
        
        indx = find(tempV == 1);
%         ts_press = leverTimeC{1,i}(indx(1));
        ts_relea = tempT(indx(1)+1);
        
        T_realign = tempT - ts_relea;
%         tempV(T_realign<0) = []; tempT(T_realign<0) = [];
        
        lever_v = tempV(T_realign<500 & T_realign > 0);
        lever_T = tempT(T_realign<500 & T_realign > 0);
        T_realign(T_realign>500 | T_realign <= 0) = [];
        lever_v(lever_v == 0) = -1;
        
        press_v(i) = sum(lever_v==1);
        release_v(i) = sum(lever_v == -1);
        
        if strcmp(input.trialOutcomeCell{1,i},'success')
            succ = succ + 1;
            press_v_s(i) = press_v(i);
            release_v_s(i) = release_v(i);
        elseif strcmp(input.trialOutcomeCell{1,i},'failure')
            fail = fail + 1;
            press_v_f(i) = press_v(i);
            release_v_f(i) = release_v(i);
        end
        
        
        
        for k = 1:length(lever_T)
            fig = plot([T_realign(k),T_realign(k)], [shift,lever_v(k)+shift]); hold on
        end
        %     if ismember(i, [5,10,15,20,25,30,35,40,45,50])
                 hline(shift);
        %     end
        shift = shift+2;
    end
end

Out.mean_p = mean(press_v);

Out.mean_r = mean(release_v);

Out.mean_p_s = sum(press_v_s)/succ;

Out.mean_r_s = sum(release_v_s)/succ;

Out.mean_p_f = sum(press_v_f)/fail;

Out.mean_r_f = sum(release_v_f)/fail;

saveas(fig, [out_dir, 'Behavior.fig']);
save([out_dir,'behavior.mat'], 'Out');
