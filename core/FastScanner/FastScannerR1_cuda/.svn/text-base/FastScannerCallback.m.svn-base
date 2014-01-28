function FastScannerCallback(obj, event)

t = event.Data.AbsTime;
fprintf(1,'%02d:%02d:%02d.%02d %s\n',t(4),t(5),floor(t(6)),...
       round((t(6)-floor(t(6)))*100),obj.Name)
