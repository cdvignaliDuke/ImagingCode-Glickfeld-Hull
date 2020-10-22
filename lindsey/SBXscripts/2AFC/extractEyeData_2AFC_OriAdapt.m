close all
clear all
clc
dataset = 'oriAdapt_V1';
eval(dataset);
nexp = size(expt,2);
LG_base = '\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\home\lindsey';

for iexp = 15:nexp
    if expt(iexp).eyeData
        mouse = expt(iexp).mouse;
        date = expt(iexp).date;
        fprintf([mouse ' ' date '\n'])
        ImgFolder = expt(iexp).runs;
        nrun = expt(iexp).nrun;
        run_str = ['runs']; 
        
        for irun = 1:nrun
            run_str = [run_str '-' expt(iexp).runs(irun,:)];
        end
        if ~exist(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
            fprintf('1')
            mkdir(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str]))
        end
        for irun = 1:nrun
            CD = [LG_base '\Data\2P_images\' expt(iexp).mouse '\' expt(iexp).date '\' expt(iexp).runs(irun,:)];
            cd(CD);
            fn = [ImgFolder(irun,:) '_000_000_eye.mat'];
            data = load(fn);          % should be a '*_eye.mat' file
            data = squeeze(data.data);      % the raw images...
            data = data(15:end-15, 40:end-40,:);
            rad_range = [10 28];
            warning off;
            A = cell(size(data,3),1);
            B = cell(size(data,3),1);
            C = cell(size(data,3),1);
            D = cell(size(data,3),1);
            for n = 1:size(data,3)
                A{n} = [0,0];
                B{n} = [0];
                C{n} = [0];
                D{n} = [0];
            end
            eye = struct('Centroid',A,'Area',B,'Val',C,'SNR',D);
            radii = [];
            Centroid = [];
            SNR = [];
            Val = [];
            Area = [];
            Eye_data = [];
            for n = 1:size(data,3)
                [center,radii,metric] = imfindcircles(squeeze(data(:,:,n)),rad_range,'Sensitivity',0.9);
                [val,idx] = max(metric);          % pick the circle with best score
                if(isempty(center))
                    eye(n).Centroid = [NaN NaN];    % could not find anything...
                    eye(n).Area = NaN;
                    eye(n).Val = NaN;
                    eye(n).SNR = NaN;
                else
                    t = double(data(:,:,n));
                    vector_of_y_values = (1:size(data,1)) - center(idx,2);
                    vector_of_x_values = (1:size(data,2)) - center(idx,1);
                    [Yg, Xg] = ndgrid(vector_of_y_values, vector_of_x_values);
                    idx1 = find(Xg.^2 + Yg.^2 < (radii(idx)/2).^2);
                    idx2 = find(Xg.^2 + Yg.^2 < (radii(idx).*2.5).^2 & Xg.^2 + Yg.^2 > (radii(idx).*1.5).^2);
                    snr = mean(t(idx1))./mean(t(idx2));
                    eye(n).SNR = snr;
                    eye(n).Val = val;
                    eye(n).Centroid = center(idx,:);
                    eye(n).Area = pi*radii(idx)^2;
                end
                if mod(n,100)==0
                    fprintf('Frame %d/%d\n',n,size(data,3));
                end
            end
            Centroid = cat(3,Centroid, cell2mat({eye.Centroid}'));
            Area = [Area cell2mat({eye.Area}')];
            Val = [Val double(cell2mat({eye.Val}'))];
            SNR = [SNR double(cell2mat({eye.SNR}'))];
            Eye_data = cat(3, Eye_data, data);
        end
        
        figure; hist(sqrt(Area./pi)); movegui('center')
        
        x1 = find(isnan(Area));
        x2 = find(~isnan(Area));
        x3 = find(Val<0.26 & SNR<1.9);

        x = unique([x1; x3]);
        if length(x)>25
            minx = 25;
        else
            minx = length(x);
        end

        frames = sort(randsample(length(x),minx));
        figure;
        start = 1;
        for i = 1:minx
            subplot(5,5,start);
            imagesq(data(:,:,x(frames(i)))); 
            hold on;
            scatter(Centroid(x(frames(i)),1), Centroid(x(frames(i)),2))
            title([num2str(chop(SNR(x(frames(i))),2)) ' ' num2str(chop(Val(x(frames(i))),2))])
            %title(num2str(x(frames(i))))
            start = start+1;
        end
        suptitle('No pupil detected')
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_noPupil.pdf']),'-dpdf','-fillpage');

        x = setdiff(x2,x3);
        frames = sort(randsample(length(x),minx));
        figure;
        start = 1;
        for i = 1:minx
            subplot(5,5,start);
            imagesq(data(:,:,x(frames(i)))); 
            hold on;
            scatter(Centroid(x(frames(i)),1), Centroid(x(frames(i)),2))
            title([num2str(chop(SNR(x(frames(i))),2)) ' ' num2str(chop(Val(x(frames(i))),2))])
            %title(num2str(x(frames(i))))
            start = start+1;
        end
        suptitle('Pupil detected')
        print(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_Pupil.pdf']),'-dpdf','-fillpage');
    end
    
    pause
    close all
    save(fullfile(LG_base, 'Analysis\2P', [date '_' mouse], [date '_' mouse '_' run_str], [date '_' mouse '_' run_str '_pupil.mat']), 'Area', 'Centroid', 'SNR', 'Val');
    
end
