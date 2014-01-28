base = 'G:\users\lindsey\analysisLG\active mice\';
date = '110518';
mouse = 'Y1';
expt = '5pos';
channel = '_green';
nsfs = 5;
ntfs = 1;
nCond = num2str(nsfs*ntfs+1);
%fn= fullfile(base, mouse, date, 'analysis', [date '_' mouse '_' expt '_reg_dF_F_Nstim' nCond '_POST_1_5.tif']);
fn = fullfile(base, mouse, date, [date '_' mouse '_' expt '.tif']);
stack = readtiff(fn);

%for 2P
[a b z] = size(stack);
stack_flat = zeros(ntfs*a, nsfs*b);

for isf = 1:nsfs
    for itf = 1:ntfs
        stack_flat(1+(itf-1)*a:a+(itf-1)*a,1+(isf-1)*b:b+(isf-1)*b) = stack(:,:,isf+(itf-1)*nsfs);
    end
end

out = fullfile(base, mouse, date, [date '_' mouse '_' expt '_SFxTF.tif']);
%out = fullfile(base, mouse, date, 'analysis', [date '_' mouse '_' expt '_ori.tif']);
writetiff(stack_flat,out);


%for old epi
[b a z] = size(stack);
stack_flat = zeros(ntfs*b, nsfs*b);

for isf = 1:nsfs
    for itf = 1:ntfs
        stack_flat(1+(itf-1)*b:b+(itf-1)*b,1+(isf-1)*b:b+(isf-1)*b) = rot90(rot90(rot90(stack(:,a-b+1:a,isf+(itf-1)*nsfs))));
    end
end

out = fullfile(base, mouse, date, [date '_' mouse '_' expt '_flip.tif']);
%out = fullfile(base, mouse, date, 'analysis', [date '_' mouse '_' expt '_ori.tif']);
writetiff(stack_flat,out);

