%vchange this to change which channel is plotted
%ch_wheel  = 1;
ch_wheel = 3;


thresh_wheelmvmt = .003;
Trigs_wheel = find(diff(data(:,ch_wheel)) > thresh_wheelmvmt);
           Trigs_wheel(find(diff(Trigs_wheel)./Fs < .01)) = [];

           SCALE_WHEEL = .048; %15 ticks per revolution, 6" diam -> .48m circumference, so .048 m/tick
            Trigs2 = Fs:Fs:size(data,1);
            Ntrigs = length(Trigs2)
            Ndown = 1;

           Wheel_mat = zeros(Ntrigs,1);
           for count = 1:Ntrigs-1
               ind = find(Trigs_wheel>(Trigs2(count)-(Fs/Ndown)) & Trigs_wheel<=(Trigs2(count)));
               %    if length(ind) > 0
               %        hi = 1
               %        pause(10);
               %    end

               Wheel_mat(count) = length(ind)*SCALE_WHEEL;
           end
           
           
           figure;
           plot(Wheel_mat);
           ylabel('wheel speed (m/s)');
           title(['Channel #',num2str(ch_wheel)]);
           xlabel('time (s)');