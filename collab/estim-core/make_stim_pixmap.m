function [outS smOutS figH] = make_stim_pixmap(stack, varargin)
% calculate dF/F pixel map
%$Id: make_stim_pixmap.m 341 2008-09-24 14:45:44Z histed $

defs = { 'StimFrNs', [], ...
         'StimFrPeakOffset', 0, ...
         'BasePreOffset', 0, ...
         'DFOFDivFrameToUse', [], ...
         'PostFrameNs', 1, ...
         'NBaseFrames', 10, ...
         'DoPlot', true, ...
         'CLim', [-0.1 0.35], ...
         'Colormap', [], ... % set below based on CLim  %cmap_posneg_yck(100, 22), ...
         'DoColorbar', true, ...
         'DoMapSmooth', true, ...
         'MapSmoothWidth', 2, ...
         'MapSmoothSigma', 3, ...
         'ComputeWhat', 'dFOF', ... % 'dF', 'stimAverage' e.g.
         'MinDenomF', 40, ...
         'DoThresholdBeforeSmooth', false, ...
         'PrintStatus', true };
     
uo = stropt2struct(stropt_defaults(defs, varargin));

% fill in cmap based on CLim
if isempty(uo.Colormap)
    % use 256-entry map
    uo.Colormap = cmap_posneg_yck(256, round(256*-uo.CLim(1)/sum(abs(uo.CLim))));    
end



%%%%

assert(uo.PostFrameNs > 0);

stC = {};
%whichStimN = 5;
%postFrNs = [10];% 2424];  % 10
%nBase = 300;
%cLim = [-0.1 0.35];
%assert(length(tsAvg.blockStimNs) == 2, 'bug');
%stimFrNs = tsAvg.stimNs(whichStimN:nStimsInSer:end);

nStims = length(uo.StimFrNs);
if uo.PrintStatus
    fprintf(1, 'Averaging pixels, %d total, stim ', nStims);
end
for iT = 1:nStims
    tStimFrN = uo.StimFrNs(iT);

    % figure out baseline frames
    assert(uo.BasePreOffset >= 0, 'use a pos val for BasePreOffset');
    startBaseFr = (tStimFrN - uo.NBaseFrames - uo.BasePreOffset);
    baseNs = startBaseFr:(startBaseFr+uo.NBaseFrames-1);
    % testing
    assert( (uo.BasePreOffset > 0  ...
             || (baseNs(end) == tStimFrN-1)), ...
            'bug in indexing: if no pre offset, last fr should be pre stimfr');
    
    % figure out stimulation frames
    stimN = (tStimFrN + uo.StimFrPeakOffset);
    stimNs = stimN:(stimN+uo.PostFrameNs-1);

    % compute average stim frame
    stimF = mean(stack(:,:,stimNs),3);
    
    switch uo.ComputeWhat
        case 'dFOF'
            baseF = mean(stack(:,:,baseNs),3);
            if isempty(uo.DFOFDivFrameToUse)
                baseD = baseF;
                if uo.MinDenomF > 0
                    baseD(baseD<uo.MinDenomF) = uo.MinDenomF;
                end
            else
                baseD = uo.DFOFDivFrameToUse;
            end
            
            
            w = warning('off');  % remove divide by zero error, user knows to deal
            stC{iT} = (stimF - baseF) ./ baseD;
            warning(w);
        case 'stimAverage'
            % just keep the average
            stC{iT} = stimF;
        case 'dF'
            baseF = mean(stack(:,:,baseNs),3);
            stC{iT} = (stimF - baseF);
        otherwise
            error('Invalid ComputeWhat: %s', uo.ComputeWhat);
    end
    
    if uo.PrintStatus
        fprintf(1,' %d', iT);
    end
end

if uo.PrintStatus
    fprintf(1, ' done.\n');
    fprintf(1, 'Computing mean and median maps... ');
end
eachDFOF = cat(3, stC{:});
%clear stC;
outS = mean(eachDFOF, 3); 
%keyboard
%outM = median(eachDFOF, 3); 
if uo.PrintStatus
    fprintf(1, 'done.\n');
end

if uo.DoMapSmooth
    if uo.DoThresholdBeforeSmooth
        topVal = uo.CLim(2);
        thrOutS = outS;
        thrOutS(thrOutS>topVal) = topVal;
    else
        thrOutS = outS;
    end
    smOutS = smooth2(thrOutS, 'gauss', ...
                     uo.MapSmoothWidth*[1 1], uo.MapSmoothSigma.*[1 1]);
else
    smOutS = outS;
end


if uo.DoPlot
    % quick map fig
    figH = figure;
    imagesc(smOutS);
    if ~isempty(uo.CLim)
        set(gca, 'CLim', uo.CLim);
    end
    colormap(uo.Colormap);
    if uo.DoColorbar
        cbH = colorbar;
        switch uo.ComputeWhat
          case 'dFOF'
            ylabel(cbH, 'dF/F');
          case 'dF'
            ylabel(cbH, 'dF');
          case stimAverage
            ylabel(cbH, 'F');
        end
    end
    set(gca, 'XTick', [], ...
             'YTick', []);
    axis square;
    title(sprintf('\\bf%s\\rm', strrep(mfilename, '_', '\_')));

else
    figH = [];
end

