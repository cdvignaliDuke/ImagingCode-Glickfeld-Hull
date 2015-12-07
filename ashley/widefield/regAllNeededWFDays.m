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
        fprintf('Starting subject %s, date %s\n', ...
        mouse, dateStr);
        eval(['i' mouse '_paths'])
        img_avg = imread(fullfile(data_folder,[dateStr '_' mouse],[mouse behav_run num2str(xd.ImgDataRun1(indexRowN))],[mouse behav_run num2str(xd.ImgDataRun1(indexRowN)) '_MMStack.ome.tif']));
        
        roi_avg = readtiff(fullfile(anal_pn, mouse, [mouse '_' roi_date '_roi_data_avg.tif']));
        AVG = double(img_avg);
        AVG(find(AVG>1e5)) = 0;
        AVG = (AVG./max(max(abs(AVG)))); % AVG can't have any depth
        target = roi_avg;
        target(find(target>1e5)) = 0;
        target = (target./max(max(abs(target)))); % Used to be plus 0.5

        [input_points, base_points] = cpselect(AVG,target,'Wait', true) 
        save(fullfile(anal_pn, mouse, [mouse '_' dateStr '_reg_data.mat']), 'img_avg', 'roi_avg', 'input_points', 'base_points')
    end

    if isempty(dateStr), break; end
end