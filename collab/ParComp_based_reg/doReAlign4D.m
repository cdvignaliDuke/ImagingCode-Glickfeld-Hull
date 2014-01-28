function RegInfo=doReAlign4D(RegInfo,trial_per_thread,FrameChunk)

datestr(clock)

si=size(RegInfo.Runs);
for i=1:si(2)
    fnums(i)=RegInfo.Runs(i).RunID;
end

Nupsamp_XYZ=RegInfo.All.Nupsamp_XYZ;
limits=RegInfo.All.limits;
nPlanes=RegInfo.All.nPlanes;
Nvolbin=RegInfo.All.Nvolbin;
is_rect=RegInfo.All.is_rect;
reverse_ch=RegInfo.All.reverse_ch;
path=RegInfo.All.path;
fbase=RegInfo.All.fbase;
ch=RegInfo.All.ch;


FrameChunk=floor(FrameChunk/nPlanes)*nPlanes;
%Align and save each run
for v=1:length(fnums)


sz=RegInfo.stack_fix_sz{v};
sz(3)=sz(3)-rem(sz(3),nPlanes);
nRead=ceil(sz(3)/FrameChunk);
fid1=fopen(fullfile(path,fbase{v},[fbase{v},ch], [fbase{v},ch,'.tif']),'r');  
fid2=fopen([path,fbase{v},'_align.bin'],'w');

for j=1:nRead
    if j~=nRead
        frames=FrameChunk;
    else
        if rem(sz(3),FrameChunk)>0
            frames=rem(sz(3),FrameChunk);
        else
            frames=FrameChunk;
        end    
        %nopar=1;
    end
    
    display(['loading original stack ',num2str(v),': up to frame ', num2str((j-1)*FrameChunk+frames),' of ', num2str(sz(3))]);
      
    stack=fread(fid1,[sz(1)*sz(2),sz(3)],'uint16=>uint16');
    stack=reshape(stack,sz(1),sz(2),sz(3));
    display('done loading original frames');
    
    nSlices = size(stack,3);
    nVol=floor(nSlices/nPlanes);
    nTrial=floor(nSlices/(Nvolbin*nPlanes));
    nChunks = ceil(nTrial/trial_per_thread);

    parameterCell = cell(1,nChunks);

    ShiftMat=RegInfo.Runs(v).ShiftMat;

    ShiftMat(find(ShiftMat(:,1) < limits(1)),1) = limits(1);
    ShiftMat(find(ShiftMat(:,1) > limits(2)),1) = limits(2);
    ShiftMat(find(ShiftMat(:,2) < limits(3)),2) = limits(3);
    ShiftMat(find(ShiftMat(:,2) > limits(4)),2) = limits(4);
    ShiftMat(find(ShiftMat(:,3) < limits(5)),3) = limits(5);
    ShiftMat(find(ShiftMat(:,3) > limits(6)),3) = limits(6);


    for i=1:nChunks
            start = (i-1)*(trial_per_thread*nPlanes*Nvolbin)+1;
            startShift = (i-1)*(trial_per_thread*Nvolbin)+1;
            if i == nChunks
                stop = nSlices;
                stopShift = nVol;
            else
                stop = start + (trial_per_thread*nPlanes*Nvolbin) - 1;
                stopShift = startShift + (trial_per_thread*Nvolbin) - 1;
            end

            %subShiftMat=ShiftMat(startShift:stopShift,:);
            subShiftMat=ShiftMat((startShift:stopShift)+((j-1)*FrameChunk/nPlanes),:);
            substack=stack(:,:,start:stop);
            parameterCell{1,i} = {substack,subShiftMat,nPlanes,is_rect,reverse_ch,Nupsamp_XYZ};
    end

    clear stack;
    
    parfor i=1:size(parameterCell,2)
         parameters=parameterCell{1,i};
         out{i} = ReAlign4D_MA3(parameters{1}, parameters{2}, parameters{3},parameters{4},...
             parameters{5},parameters{6});
    end  
    

    for i=1:(nChunks-1)
        stack4D_reg(:,:,:,:,i)=out{1,i};
    end

    si=size(stack4D_reg);
    stack4D_reg=reshape(stack4D_reg,si(1),si(2),si(3),si(4)*si(5));
    si=size(stack4D_reg);
    si2=size(out{1,nChunks});
    stack4D_reg(:,:,:,(si(4)+1):(si(4)+si2(4)))=out{1,nChunks};
    
    disp('Subsection alignment complete.');
    
    clear out

    stack4D_reg_avg_up = squeeze(mean(stack4D_reg,4));
    stack4D_vec_avg_up_sub(:,:,:,j)=stack4D_reg_avg_up;

    stack4D_reg_max_up = squeeze(max(stack4D_reg,[],4));
    stack4D_vec_max_up_sub(:,:,:,j)=stack4D_reg_max_up;
    
    %downsample back to original frame #
    stack4D_reg=stack4D_reg(:,:,[1:ceil(si(3)./2)]*2-1,:);
    volsz=size(stack4D_reg);
    
    %write to file
    fwrite(fid2,uint16(stack4D_reg),'uint16');

    clear stack4D_reg
    clear stack4D_reg_avg_up
    clear stack4D_reg_max_up
    
    disp('Subsection output complete.');
end

fclose(fid2);
fclose(fid1);
RegInfo.stack_align_sz{i}=volsz(1:3);

stack4D_reg_avg_up = squeeze(mean(stack4D_vec_avg_up_sub,4));
b=writetiff(stack4D_reg_avg_up,[path,fbase{v},'_align_avg_up.tif']);
stack4D_vec_avg_up(:,:,:,v)=stack4D_reg_avg_up;

stack4D_reg_avg=stack4D_reg_avg_up(:,:,[1:ceil(si(3)./2)]*2-1);
b=writetiff(stack4D_reg_avg,[path,fbase{v},'_align_avg.tif']);
stack4D_vec_avg(:,:,:,v)=stack4D_reg_avg;

stack4D_reg_max_up = squeeze(max(stack4D_vec_max_up_sub,[],4));
b=writetiff(stack4D_reg_max_up,[path,fbase{v},'_align_max_up.tif']);
stack4D_vec_max_up(:,:,:,v)=stack4D_reg_max_up;

stack4D_reg_max=stack4D_reg_max_up(:,:,[1:ceil(si(3)./2)]*2-1);
b=writetiff(stack4D_reg_max,[path,fbase{v},'_align_max.tif']);
stack4D_vec_max(:,:,:,v)=stack4D_reg_max;

clear stack4D_reg
clear stack4D_reg_avg_up
clear stack4D_reg_avg
clear stack4D_reg_max_up
clear stack4D_reg_max

disp(['Run',num2str(fnums(v)),' output complete.']);
end

RegInfo.All.avg_up=squeeze(mean(stack4D_vec_avg_up,4));
RegInfo.All.avg=squeeze(mean(stack4D_vec_avg,4));
RegInfo.All.max_up=squeeze(mean(stack4D_vec_max_up,4));
RegInfo.All.max=squeeze(mean(stack4D_vec_max,4));

alloutfname=RegInfo.All.alloutfname;
b=writetiff(RegInfo.All.avg_up,[alloutfname,'align_avg_up.tif']);
b=writetiff(RegInfo.All.avg,[alloutfname,'align_avg.tif']);
b=writetiff(RegInfo.All.max_up,[alloutfname,'align_max_up.tif']);
b=writetiff(RegInfo.All.max,[alloutfname,'align_max.tif']);

RegInfo.All.limits=limits;
RegInfo.All.Nupsamp_XYZ=Nupsamp_XYZ;
RegInfo.All.fnums=fnums;

save([alloutfname,'RegInfo.mat'],'RegInfo');

datestr(clock)