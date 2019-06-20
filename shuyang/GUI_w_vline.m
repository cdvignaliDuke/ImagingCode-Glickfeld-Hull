%used in twoP_MDspikes. same thing as function GUI. should be used after spike threshold is determined.
%This function has vertical lines at certain time points to represent
%onset/offset, everything else is the same as GUI.
%generate a GUI with plots of raw fluorescnece for cells during picked
%frames (this can be run triggered ave, etc). 
%Spike events (Fluorescnece > std best) have red dots over them.
function [fig] = GUI(TCave, frames_mat,spk_inx)
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
        r = randperm(size(frames_mat,2),4); 
        for j = 1:4                                                       % first 9 running trig windows
            frm_abv = intersect(spk_inx{k},frames_mat(:,r(j))); %get the frames when spike happens
            subplot(2,2,j);plot(TCave(frames_mat(:,r(j)),k));hold on; vline(31,'r'); hold on;
            plot(frm_abv - min(frames_mat(:,r(j))), (max(TCave(frames_mat(:,r(j)),k))-10)*ones(1,length(frm_abv)),'r.'); % plot red dots on top of the spikes
            % do frm_abv - min because otherwise the x of the red dots start from 0.
            set(gca,'xticklabel',[]); %,'yticklabel',[]);
            title(['cell' num2str(k) 'frm' num2str(frames_mat(1,r(j))) '-' num2str(frames_mat(45,r(j)))]);
            drawnow; hold off
        end
        
    end
    
end


