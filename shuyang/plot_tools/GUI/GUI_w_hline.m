% used in twoP_MDspikes, helpful for determining spike threshold.
% should be used after spike threshold is determined.
% generate a GUI with plots of raw fluorescence for each cell during
% several time periods (these time periods are fixed, can't drag and see other times ). 
% hlines represent the std, 2std, 2.5std, 3std of derivatives.
function [fig] = GUI_w_hline(deriv, frames_mat,std_deriv,std2,std2_5,std3)
num_components = size(deriv,2);
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
        for j = 1:4                                                         % every time period
            subplot(2,2,j);plot(deriv(frames_mat(:,j),k));hold on; 
            hline(std_deriv(k),'r'); hold on;
            hline(std2(k),'b');hold on;
            hline(std2_5(k),'g');hold on;
            hline(std3(k),'k');hold on;
            set(gca,'xticklabel',[]); %,'yticklabel',[]);
            title(['cell' num2str(k) 'frm' num2str(frames_mat(1,j)) '-' num2str(frames_mat(end,j))]);
            drawnow; hold off
        end
        
    end
    
end


