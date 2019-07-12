%used in twoP_MDspikes. should be used after spike threshold is determined.
%generate a GUI with plots of raw fluorescnece for cells during picked
%frames (this can be run triggered ave, etc). 
%Spike events (Fluorescnece > std best) have red dots over them.
function [fig] = GUI(TCave, frames_mat,spk_inx,sessions)
num_components = size(TCave,2);
ind_all = 1:size(TCave,1);
fig = figure('Visible','off');
set(gcf,'Position',2*[300 300 960 480]);
set(gcf,'PaperPosition',2*[300 300 960 480]);

sld = uicontrol('Style','slider',...
    'Min',1,'Max',num_components,'Value',1,'SliderStep',[1/(num_components-1) 1],...
    'Position',[150 20 800 20],...
    'Callback',@surfzlim); %num_components is number of total images

txt = uicontrol('Style','text',...
    'Position',[400 45 120 20],...
    'String','Component');

fig.Visible = 'on';
plot_component(1)

    function surfzlim(source, callbackdata)
        k = source.Value;
        plot_component(round(k))
    end

    function plot_component(k)                                              % each cell
    cla
    %for i = randperm(size(TCave,2),10)                                      
        %figure('units', 'normalized', 'outerposition', [0 0 1 1]); % open the figure full screen
        %set(gcf, 'position', get(0,'screensize'))
        for j = 1:4                                                         % first 4 windows
            %subplot(2,2,j);plot(min(frames_mat(:,j)):max(frames_mat(:,j)),TCave(frames_mat(:,j),k));
            subplot(2,2,j);plot(min(frames_mat(:,j)):max(frames_mat(:,j)),TCave(frames_mat(:,j),k));

            frm_abv = intersect(spk_inx{k},frames_mat(:,j)); %get the frames when spike happens
            frames = frames_mat(:,j);
            spklogic = ismember(frames,frm_abv);
            plotred = TCave(frames_mat(:,j),k).*spklogic;
            plotred(plotred==0) = NaN;
            plotx = frames_mat(:,j).*spklogic;
            plotx(plotx ==0) = NaN;
            hold on;
            plot(plotx,plotred,'ro','MarkerSize',8);
            %plot(frm_abv, (max(TCave(frames_mat(:,j),k))-10)*ones(1,length(frm_abv)),'r.'); % plot red dots on top of the spikes
            xlim([min(frames_mat(:,j)),max(frames_mat(:,j))]);
            hold on; set(gca,'xticklabel',[]); %,'yticklabel',[]);
            hold on; title(['cell' num2str(k) 'frm' num2str(frames_mat(1,j)) '-' num2str(frames_mat(end,j))]);
            drawnow; hold off
        end
    end
        supertitle(['TC ave_' sessions]); hold off
   
    
end


