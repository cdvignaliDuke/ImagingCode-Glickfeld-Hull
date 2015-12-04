function regAllNeededWFDays(mouse)
while true
    close all;
    rc = behavConstsWF(mouse);
    
    [xd indexRowN] = findNextIndexForRegWF(rc);
    
    if isempty(indexRowN)
        disp('No more records to do');
        dateStr = [];
        return
    else
        % do the one we found
        dateStr = xd.DateStr{indexRowN};
        fprintf('Starting subject %s, date %s', ...
        mouse, dateStr);
        eval(['i' mouse '_paths'])
        pn_img_folder = dir(fullfile(data_folder,[dateStr '_' mouse],[mouse behav_run '*']));
        pn_img_data = dir(fullfile(data_folder,[dateStr '_' mouse],pn_img_folder.name(1,:),[mouse behav_run '*']));
        %img_data = readtiff(fullfile(data_folder,[dateStr '_' mouse],pn_img_folder.name(irun,:), pn_img_data.name(irun,:)));

        %img_avg = mean(img_data, 3);
        load('Z:\home\lindsey\Analysis\Widefield_imaging\633\633_151103_reg_data.mat')
        roi_avg = readtiff(fullfile(anal_pn, mouse, [mouse '_' roi_date '_roi_data_avg.tif']));
        AVG = img_avg;
        AVG(find(AVG>1e5)) = 0;
        AVG = (AVG./max(max(abs(AVG)))); % AVG can't have any depth
        target = roi_avg;
        target(find(target>1e5)) = 0;
        target = (target./max(max(abs(target))))-.5; % Used to be plus 0.5

        [input_points, base_points] = cpselect(AVG,target,'Wait', true) 
        save(fullfile(anal_pn, mouse, [mouse '_' dateStr '_reg_data.mat']), 'img_avg', 'roi_avg', 'input_points', 'base_points')
    end

    if isempty(dateStr), break; end
end