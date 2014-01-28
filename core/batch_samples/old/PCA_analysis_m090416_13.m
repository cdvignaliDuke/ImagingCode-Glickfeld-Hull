fname='E:\users\kenichi\leicatemp\rtifpc\rvstim13pc.tif'
outfname='090416_13';
N_PC=10; % Number of PCs to extract
Noff=10;
Non=10;
Nstim=12;
ker=fspecial('gaussian',10,1.5);


%%%%%%%%%%%%%%
% PCA decomposition

Nframe_per_run=(Non+Noff)*Nstim;
stack=readtiff(fname);
[Nx,Ny,Nt]=size(stack);
Nrep=floor(Nt/Nframe_per_run);
Nframe=Nframe_per_run*Nrep;
stack2=double(reshape(stack,Nx*Ny,Nt));
stack2=stack2(:,1:Nframe);
avg=squeeze(mean(stack2,2));
stack2=stack2-repmat(avg,1,Nframe);
[U,S,V] = pca(stack2,N_PC);
U2=reshape(U,Nx,Ny,N_PC);
save PCAresult U S V U2

%%%%%%%%%%%%%%%
% display PCA timecourses

figure;
for i=1:N_PC
    subplot(N_PC,3,i*3-2);
    plot(V(:,i));
    subplot(N_PC,3,i*3-1);
    plot(sum(reshape(V(:,i),Nframe_per_run,Nrep),2));
    set(gca,'xtick',[Noff:Noff+Non:Nframe_per_run-Non]);
    set(gca,'xgrid','on');
    subplot(N_PC,3,i*3);
    power=abs(fft(V(:,i)));
    plot([0;power(2:ceil(Nframe/2))]);
end
saveas (gcf, 'PCA_tcourses.fig');

%%%%%%%%%%%%%%%
% display PCA images

figure;
for i=1:N_PC
    subplot(4,3,i);
    imagesc(U2(:,:,i));
    colormap(gray);
end

saveas (gcf, 'PCA_images.fig');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   noise removal
%

nstack=stack2+repmat(avg,1,Nframe);

    nstack2=reshape(nstack,Nx,Ny,Nframe);
    RMVname=[outfname, '_RMV0'];
    nstack_avg=mean(reshape(nstack2,Nx,Ny,Nframe_per_run,Nrep),4);
    impixtimecourse3(nstack_avg,4,1,[0.95,1.1],[Noff:Noff+Non:Nframe_per_run-Non],RMVname);
    save([RMVname,'.mat'],'nstack_avg');
%      nstack_avg2=mean(reshape(nstack_avg(:,:,[240,1:239]),512,512,12,20),3);
%      base=mean(nstack_avg(:,:,[11:24:11+24*9]),3);
%      dF=nstack_avg2(:,:,[2:2:20])-repmat(base,[1 1 10]);
%      for i=1:10
%          subplot(3,4,i);
%          imshow(filter2(ker,dF(:,:,i))/100);
%      end
%      saveas (gcf, [RMVname,'_dFimages.fig']);


RMV=N_PC;
for i=1:RMV
    nstack=nstack-U(:,i)*S(i,i)*(V(:,i)-repmat(mean(V(:,i),1),Nframe,1))';
    nstack2=reshape(nstack,Nx,Ny,Nframe);
    RMVname=[outfname, '_RMV', num2str(i)];
    writetiff(nstack2,[RMVname,'.tif']);
    nstack_avg=mean(reshape(nstack2,Nx,Ny,Nframe_per_run,Nrep),4);
    figure;
    impixtimecourse3(nstack_avg,4,1,[0.95,1.1],[Noff:Noff+Non:Nframe_per_run-Non],RMVname);
    save([RMVname,'.mat'],'nstack_avg');
%      nstack_avg2=mean(reshape(nstack_avg(:,:,[240,1:239]),512,512,12,20),3);
%      base=mean(nstack_avg(:,:,[11:24:11+24*9]),3);
%      dF=nstack_avg2(:,:,[2:2:20])-repmat(base,[1 1 10]);
%      for i=1:10
%          subplot(3,4,i);
%          imshow(filter2(ker,dF(:,:,i))/100);
%      end
%      saveas (gcf, [RMVname,'_dFimages.fig']);
end