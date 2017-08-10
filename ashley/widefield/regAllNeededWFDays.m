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
        data_folder = fullfile(rc.dataRootDir, [dateStr '_' mouse]);

        img_avg = imread(fullfile(data_folder,[dateStr '_' mouse '_1'], [dateStr '_' mouse '_1_MMStack.ome.tif']));
        
        if indexRowN == 1
            roi_avg = img_avg;
            writetiff(roi_avg, fullfile(rc.structOutput,[mouse '_roi_avg.tif']))
            input_points = ones(3,2);
            base_points = ones(3,2);
        else
            roi_avg = readtiff(fullfile(rc.structOutput,[mouse '_roi_avg.tif']));
            AVG = double(img_avg);
            AVG(find(AVG>1e5)) = 0;
            AVG = (AVG./max(max(abs(AVG)))); % AVG can't have any depth
            target = double(roi_avg);
            target(find(target>1e5)) = 0;
            target = (target./max(max(abs(target)))); % Used to be plus 0.5

            [input_points, base_points] = cpselect(AVG,target,'Wait', true) 
        end
        save(fullfile(rc.structOutput, [mouse '_' dateStr '_reg_data.mat']), 'img_avg', 'roi_avg', 'input_points', 'base_points')
    end

    if isempty(dateStr), break; end
end