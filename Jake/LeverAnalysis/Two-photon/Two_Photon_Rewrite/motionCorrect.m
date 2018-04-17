function [img_reg, laser_on_ind_conserv, nt] = motionCorrect(data_dir, runID, out_dir)
[img, ~, ~] = loadFile( data_dir, runID, [], []);
nt = size(img,3);
useGPU = 0;
usFacs = 10;
% [npw, nph, nframes1] = size(img);

if exist([out_dir, 'Reg_out.mat'],'file') == 2
    out1 = load([out_dir, 'Reg_out.mat']);
    img_reg = [];
    if useGPU == 1
        for iF = 1 : floor(nframes1/1000)
            [~,img_reg_temp]=stackRegister_MA(gpuArray(img(:,:,1 : 1000)),gpuArray(img(:,:,1)),usFacs,...
                gpuArray(out1.reg_out(1 : 1000,:)));
            img(:,:,1:1000) = [];
            img_reg = cat(3, img_reg, gather(img_reg_temp));
            clear img_reg_temp
        end
        %                 elseif useWarp == 1
        %                     img_reg = zeros(size(img));
        %                     for iF = 1:nframes1
        %
        %                         [~, img_reg(:,:,iF)] = imregdemons(img(:,:,1),out1.img_ref, [32 16 8 4],'AccumulatedFieldSmoothing',2.5,'PyramidLevels',4);
        %                         img(:,:,1) = [];
        %                     end
        
    else
        [~,img_reg]=stackRegister_MA(img,img(:,:,1),usFacs,double(out1.reg_out));
        clear img
    end
else
    mkdir(out_dir)
    config_fn = dir(fullfile(data_dir,['*' runID '.mat']));
    laser_power_fn = dir(fullfile(data_dir,['*' runID '_realtime.mat']));
    [img_mat_file, laser_power_vec_ttl] = get_laser_power_data(data_dir, config_fn, laser_power_fn);
    if isempty(laser_power_vec_ttl)
        laser_power_vec_ttl = ones(1,nt);
        frame_nums_for_ref30 = randi([1,nt],1,30);
        
        sf = floor((nt - 11)/100);
        frame_nums_for_samp100 = 11:sf:nt;
        laser_on_ind_conserv = 1:nt;
    else
        laser_on_ind_conserv = conservative_laser_on(laser_power_vec_ttl);
        frame_nums_for_ref30 = laser_on_ind_conserv(randi([1,length(laser_on_ind_conserv)],1,30));
        frame_nums_for_samp100 = laser_on_ind_conserv(round(linspace(1,length(laser_on_ind_conserv))));
    end
    laser_on_ind = find(laser_power_vec_ttl);   %frame numbers of frames with laser power on determined by _realtime ttl_log
    ref30 = img(:,:,frame_nums_for_ref30);
    samp100 = img(:,:,frame_nums_for_samp100);
    dshift = [];
    for r = 1:size(ref30,3)
        [reg_out,aa] = stackRegister(samp100, ref30(:,:,r));
        dshift = [dshift;mean(((reg_out(:,3).^2)+(reg_out(:,4).^2)))];
    end
    
    %pick the frame which had the lowest dshift and motion register full movie to that frame
    min_f = find(dshift == min(dshift));
    img_ref = ref30(:,:,min_f);
    [reg_out, img_reg] = stackRegister(img(:,:,laser_on_ind_conserv), img_ref);
%     img_reg = img_reg(:,:,laser_on_ind_conserv);
    save([out_dir, 'img_reg.mat'], 'img_reg', '-v7.3');
    save([out_dir, 'Reg_out.mat'], 'reg_out','img_ref', 'laser_on_ind_conserv', 'nt');
    
    
end