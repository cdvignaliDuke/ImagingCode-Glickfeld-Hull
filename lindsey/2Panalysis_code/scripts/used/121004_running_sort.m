clear all
areas = ['PM'; 'LM'; 'Al'; 'RL'; 'AM'];
inj = 'V1';
P = 2;
Fs_fastrig = 30.82; %Hz
dTrig0 = 1/Fs_fastrig;
Nframes_pre = 144;
Nframes_post = 144;
ch_wheel = 1;

for iArea = 1:5
    matrix = 'SF5xTF5';
    image = areas(iArea,:);
    sum_base = 'G:\users\lindsey\analysisLG\experiments';
    list_fn = fullfile(sum_base, [matrix '_' num2str(P) 'P_' image], [matrix '_' num2str(P) 'P_' inj '_' image '_exp_list.mat']);
    load(list_fn);
    nexp = size(exp_list.mouse_mat,2);    
    for iexp = 1:nexp
        mouse = char(exp_list.mouse_mat{iexp});
        date = char(exp_list.date_mat{iexp});
        userun = exp_list.runs_mat{iexp};
        count_protocol = exp_list.prot_mat{iexp};
        run = exp_list.run_mat{iexp};
        blanks = exp_list.blanks_mat{iexp};
        dirs = exp_list.dir_mat{iexp};
        base = 'G:\users\lindsey\analysisLG\active mice';
        outDir = fullfile(base, mouse, date);

        if run == 1
            if dirs ==1
                nCond = 25;
            elseif dirs ==2
                nCond = 50;
            end
            DIR_base = fullfile('\\zoloto\bigstorlab\fastrig\running\Running_data',[date '_' mouse '\']);
            eval(['PARAMS_' date '_' mouse]);
            data_base = 'G:\users\lindsey\analysisLG\active mice';
            Wheel_mat = [];
            for iRun  = 1:length(userun);
                %figure out #frames
                file_USE = deblank(file_mat(userun(iRun),:));
                Nframes1 = sizetiff(fullfile(data_base,mouse,date,file_USE))*12;               

                %load running data
                DIR_data = dir([DIR_base, '*' [num2str(userun(iRun)) '.mat']]);
                load(fullfile(DIR_base,DIR_data.name));
                Ntimestamps0 = size(data,1);

                %ID trigs
                Trigs = ceil(Fs*[dTrig0:dTrig0:(dTrig0*Nframes1)]);

                %Find eye mvmt and wheel movement in the same period
                thresh_wheelmvmt = 3;
                Trigs_wheel = find(diff(data(:,ch_wheel)) > thresh_wheelmvmt);
                Trigs_wheel(find(diff(Trigs_wheel)./Fs < .01)) = [];

                SCALE_WHEEL = .048; %10 ticks per revolution, 6" diam -> .48m circumference, so .048 m/tick

                Nstim = Nframes1./(Nframes_pre+Nframes_post);
                Wheel_mat_temp = zeros(Nstim,1);
                for count = 1:Nstim-1
                   ind = find(Trigs_wheel>(Trigs(1+Nframes_pre+((count-1)*(Nframes_pre+Nframes_post)))) & Trigs_wheel<=(Trigs(count*(Nframes_pre+Nframes_post))));
                   Wheel_mat_temp(count) = length(ind)./(Nframes_post/Fs_fastrig) * SCALE_WHEEL;
                end

                Wheel_mat = [Wheel_mat; Wheel_mat_temp];
            end
            clear('run_data');

            %resort wheel data
            seqfile = [date '_' mouse '_run' num2str(userun) '_Big_Seqposition.mat'];
            load(fullfile(outDir,'analysis',seqfile));

            Wheel_sorted = zeros(size(Wheel_mat));
            start = 1;
            for iCond = 1:nCond;
                nRep = length(Big_Seqposition(iCond).ind);
                for iRep = 1:nRep;
                    ind = Big_Seqposition(iCond).ind(iRep);
                    Wheel_sorted(start) = Wheel_mat(ind);
                    start = start+1;
                end
            end

            nblanks = length(Big_Seqposition(end).ind);
            for iblank = 1:nblanks;
                ind = Big_Seqposition(end).ind(iblank);
                Wheel_sorted(start) = Wheel_mat(ind);
                start = start+1;
            end
            fn_out = fullfile(outDir,'analysis', [date '_' mouse '_run' num2str(userun) '_wheel.mat']);
            save(fn_out, 'Wheel_mat', 'Wheel_sorted');
        end
    end
end