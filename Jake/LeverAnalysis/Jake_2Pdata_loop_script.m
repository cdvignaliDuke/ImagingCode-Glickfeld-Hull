clear all
base_dir =  '\\CRASH.dhe.duke.edu\data\home\jake\Data\2P_imaging\';
date_mat = strvcat('150703', '150703','150704','150704');
run_mat = strvcat('_000_000', '_000_001', '_000_001', '_000_000');
mouse_mat = strvcat('img24','img25','img24','img25');
subNum_mat = strvcat('924','925','924','925');
time_mat = strvcat('1821', '2005','1845','1937');
out_base = 'Z:\home\lindsey\Analysis\2P\Jake\';
tc_type = 'ICAsig';
for id = 1:size(date_mat,1)
    date = date_mat(id,:);
    run = run_mat(id,:);
    mouse = mouse_mat(id,:);
    subNum = subNum_mat(id,:);
    time = time_mat(id,:);
    
    run_name = [date '_' mouse '_run' run(length(run)-2:end)];
    out_path = fullfile(out_base,run_name);
    dest =  fullfile(out_path,run_name);

    img = readtiff([dest '_ROI.tif']);

    %prep for pca
    global stack
    stack = single(img);
    defaultopts = {'nComp',300,'BorderWidth',4};
    options = cell2struct(defaultopts(2:2:end),defaultopts(1:2:end),2);
    [ny,nx,nt]=size(stack);
    roi = imborder([ny,nx],options.BorderWidth,0); 
    fprintf('Masking edges... ');
    stack= bsxfun(@times,stack,single(roi));
    fprintf('Done\n');
    % compute thin svd using randomized algorithm
    pcs = stackFastPCA(1,options.nComp);

    PCuse = [1:100];
    mu = 0;
    nIC = 32;
    termtol = 1e-6;
    maxrounds = 400;
    mixedsig = pcs.v';
    mixedfilters = pcs.U;
    CovEvals = diag(pcs.s).^2;
    [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, ...
        mixedfilters, CovEvals, PCuse, mu, nIC,[],termtol,maxrounds);
    
    save([dest '_ICs_fullres.mat'],'sm', 'ica_sig');
    
    data_tc = ica_sig';
    HAD_cmp_success_fail_2P_ROIs
end
