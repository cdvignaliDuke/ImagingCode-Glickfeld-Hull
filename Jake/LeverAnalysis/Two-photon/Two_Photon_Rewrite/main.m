% main
clear
file_info;
usFac = 4;
for i = 1:size(mouseID,2)
    for runID = 1:2
        out_dir  = fullfile('C:','Users','ziye','Documents','MATLAB','2P_Analysis',mouseID{i});
        [img_raw, skip_run] = loadFile(i);
        
        if skip_run == 1
            break
        end
        nFrames = size(img_raw,3);
        img_reg = 0.*(img_raw);
        fm = 1;
        % motion correction 
        while  fm < nFrames
            if mod(fm,1000) == 0
                fprintf('registering %i frames...\n', fm);
            end
            if fm == 1
                fprintf('registering second frame...\n');
                img_reg(:,:,fm) = img_raw(:,:,fm);
                [outs,img_reg(:,:,fm+1)]=stackRegister(img_raw(:,:,fm+1),img_reg(:,:,fm),usFac);
            else
                [outs,img_reg(:,:,fm+1)]=stackRegister(img_raw(:,:,fm+1),img_reg(:,:,fm),usFac);
            end
            fm = fm + 1;
        end
        % smooth
        img_sm = stackFilter(img_reg);
        
        
    end
end
