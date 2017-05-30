function ICbad = ChooseICs(icasig)
%ind = 1;
n = 10;
nmask = size(icasig,3);
ICbad = [];
for k  = 1:nmask/n
    for ic = (k-1)*n+1 : k*n
        subplot(n/2,2,ic);                 %change here too
        %imstretch(sm(:,:,ic),[.5 .99],1.5);
        imagesc(icasig(:,:,ic));
        %ind = ind+1;
        text(.8,.1,num2str(ic),'fontsize',12,'color','w','fontweight','bold','unit','norm');
        ICbad_input = input('Number of bad IC ', 's');
        ICbad = [ICbad str2num(ICbad_input)];
    end
end