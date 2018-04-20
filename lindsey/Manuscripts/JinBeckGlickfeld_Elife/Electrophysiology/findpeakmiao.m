function [peak,latency,peakIX]=findpeakmiao(meanTrace,timewindow, on_edges,basebins)



[M,I] = max(meanTrace(timewindow));
if  M~=0
    peak=M;
    latency= on_edges(I)*1000;
    peakIX= I+basebins; 
else
   peak=NaN;
   latency= NaN;
   peakIX= NaN; 
end


end