%used in spike_derivative_threshold. should be used after spike threshold is determined.
%generate a GUI with subtplots:
% 1. plots of raw fluorescnece for cells. Spike events (Fluorescnece > std best) have red dots over them.
% 2. plots of first derivatives of all cells

function [fig] = GUI_rawTrace_nderivThresh(TCave,deriv,spk_bi_cellmat,sessions)
std_deriv = std(deriv);
std2 = 2*std_deriv;
std3 = 3*std_deriv;
std2_5 = 2.5*std_deriv;
num_components = size(TCave,2);
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
        plotred = TCave(:,k).*spk_bi_cellmat(:,k);
        plotred(plotred==0) = NaN;
        subplot(2,1,1)
        plot(TCave(:,k)); hold on;
        plot(plotred,'ro','MarkerSize',8); ylabel('TCave');
        xlim([1,1500]);
        lims = ylim;
        text(10,lims(2)-20,['cell' num2str(k)]);
        
        drawnow; hold off
        
        subplot(2,1,2)
        plot(deriv(:,k)); hold on;
        hline(std_deriv(k),'r'); hold on;
        hline(std2(k),'b');hold on;
        hline(std2_5(k),'g');hold on;
        hline(std3(k),'k');hold on;
        xlim([1,1500]);ylabel('deconvolved signal');xlabel('frame');
        drawnow; hold off
        
    end

supertitle(['TC ave and derivative ' sessions]); hold off


end


