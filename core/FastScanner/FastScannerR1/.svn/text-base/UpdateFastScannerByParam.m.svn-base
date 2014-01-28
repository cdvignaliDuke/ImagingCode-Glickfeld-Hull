function UpdateFastScannerByParam;

global fs;

%fs.pCell.max
set(fs.handles.sliderPCellMax,'Value',fs.pCell.max);
FastScanner('sliderPCellMax_Callback',fs.handles.sliderPCellMax, [], fs.handles);

%fs.pCell.max
set(fs.handles.sliderPCellMin,'Value',fs.pCell.min);
FastScanner('sliderPCellMin_Callback',fs.handles.sliderPCellMin, [], fs.handles);

if fs.stage.JS==0
    %fs.position.x
    set(fs.handles.sliderXPos,'Value',fs.position.x);
    FastScanner('sliderXPos_Callback',fs.handles.sliderXPos, [], fs.handles);

    %fs.position.y
    set(fs.handles.sliderYPos,'Value',fs.position.y);
    FastScanner('sliderYPos_Callback',fs.handles.sliderYPos, [], fs.handles);

    %fs.position.z
    set(fs.handles.sliderZPos,'Value',fs.position.z);
    FastScanner('sliderZPos_Callback',fs.handles.sliderZPos, [], fs.handles);
end

%fs.position.zoom
set(fs.handles.sliderZoom,'Value',fs.position.zoom);
FastScanner('sliderZoom_Callback',fs.handles.sliderZoom, [], fs.handles);

%fs.cycle.XStart
set(fs.handles.txtStartX,'String',fs.cycle.XStart);

%fs.cycle.YStart
set(fs.handles.txtStartY,'String',fs.cycle.YStart);

% fs.cycle.ZStart
set(fs.handles.txtStartZ,'String',fs.cycle.ZStart);

% fs.cycle.XStep
set(fs.handles.txtStepX,'String',fs.cycle.XStep);

% fs.cycle.YStep
set(fs.handles.txtStepY,'String',fs.cycle.YStep);

% fs.cycle.ZStep
set(fs.handles.txtStepZ,'String',fs.cycle.ZStep);

% fs.cycle.NSteps
set(fs.handles.txtNofSteps,'String',fs.cycle.NSteps);

% fs.cycle.NReps
set(fs.handles.txtNofRepetitions,'String',fs.cycle.NReps);

% fs.cycle.NFramesPrStep
set(fs.handles.txtFramesPerStep,'String',fs.cycle.NFramesPrStep);

% xxdisabled: fs.cycle.Avg
set(fs.handles.chkAvg,'Value',fs.cycle.Avg);
FastScanner('chkAvg_Callback',fs.handles.chkAvg, [], fs.handles);
set(fs.handles.chkAvg, 'Enable', 'off'); % disabled

% fs.pCell.LUTfile
set(fs.handles.txtPCellLUTFile,'String',fs.pCell.LUTfile);


set(fs.handles.txtBaseFileName,'String',fs.DAQ.BaseFileName);

if isempty(fs.pCell.LUT)
    set(fs.handles.lblpCellUnits,'String','V');
    set(fs.handles.lblpCellUnits2,'String','V');
else
    set(fs.handles.lblpCellUnits,'String','mW');
    [cmax, indmax]=max(fs.pCell.LUT(:,2));
    [cmin, indmin]=min(fs.pCell.LUT(:,2));
    set(fs.handles.sliderPCellMax,'Min',cmin);
    set(fs.handles.sliderPCellMax,'Max',cmax);

    set(fs.handles.lblpCellUnits2,'String','mW');
    set(fs.handles.sliderPCellMin,'Min',cmin);
    set(fs.handles.sliderPCellMin,'Max',cmax);

    if fs.pCell.max<cmin
        fs.pCell.max=cmin;
    end
    if fs.pCell.max>cmax
        fs.pCell.max=cmax;
    end
    if fs.pCell.min<cmin
        fs.pCell.min=cmin;
    end
    if fs.pCell.min>cmax
        fs.pCell.min=cmin;
    end

    set(fs.handles.sliderPCellMax,'Value',fs.pCell.max);
    set(fs.handles.txtPCellMax,'String',num2str(fs.pCell.max));
    set(fs.handles.sliderPCellMin,'Value',fs.pCell.min);
    set(fs.handles.txtPCellMin,'String',num2str(fs.pCell.min));

    set(fs.handles.sliderPCellMin,'SliderStep',[1/(cmax-cmin) 10/(cmax-cmin)]);
    set(fs.handles.sliderPCellMax,'SliderStep',[1/(cmax-cmin) 10/(cmax-cmin)]);
end

return;
