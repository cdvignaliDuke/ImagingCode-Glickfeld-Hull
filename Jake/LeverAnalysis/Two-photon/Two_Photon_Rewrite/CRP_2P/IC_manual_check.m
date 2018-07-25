function IC_use = IC_manual_check(icasig);
%function with user input for selecting which ICs to use based upon the
%morphology in their IC intensity masks.

%assign output variable and create figure 
IC_use = [];
figure('rend', 'painters', 'pos', [50 150 1200 475]);

for ICnum = 1:size(icasig,3);
    %plot this IC
    imagesc(icasig(:,:,ICnum));
    title(['IC # ', num2str(ICnum),' of ', num2str(size(icasig,3)), ':', 'Press "y" to keep, "n" to reject.']);
    
    %get user input on this IC and log decision in output variable.
    exit_cond = 0;
    while exit_cond == 0;
        user_input = input('', 's');
        if strcmp(user_input, 'y');
            IC_use(ICnum) = 1;
            exit_cond = 1;
        elseif strcmp(user_input, 'n');
            IC_use(ICnum) = 0;
            exit_cond = 1;
        else
            disp('Invalid user input. Please troubleshoot user and try again.')
        end
    end
end
disp('IC selection complete');
close 
return
