clear all; close all; clc 

%% set parameters

%This file requires binfileload.m, autospec.m, specgram.m,
%UW_Sensitivity.m, and the 4 hydrophone sensitivity csv files in order to
%work.

%Choose whether you want to analyze all files from the folder or only a
%select few. Set all_files to true to analyze all files. If set to false,
%only file numbers in IDnum will be analyzed.
all_files = false;
num_files = 20;
IDnum = [2, 3, 4, 11];

% If all_files is false and the files selected are to be cross-correlated
% and time averaged, XandA should be set to true.
XandA = false;

% recording parameters
path = 'D:\2020-2-20 Tank range to wall test';
t_rec = 5.5;
fs = 1240000;
N = fs*t_rec;
ns = 2^15;

% this script assumes an input channel in AI0 and an output channel in AI1

%% read in data

if all_files == true
    input = zeros(num_files, N);
    output = zeros(num_files, N);
    
    for n = 1:num_files
        input(n,:) = binfileload(path,'ID', n, 0);
        output(n,:) = binfileload(path,'ID', n, 1);

    end
    
else
    input = zeros(length(IDnum), N);
    output = zeros(length(IDnum), N);

    for n = 1:length(IDnum)
        input(n,:) = binfileload(path,'ID', IDnum(n), 0);
        output(n,:) = binfileload(path,'ID', IDnum(n), 1);

    end
end

disp('Data Loaded')


%% autospectra
if all_files == true
    Gin = zeros(floor(N/ns), floor(ns/4), num_files);
    Gout = zeros(floor(N/ns), floor(ns/4), num_files);

    for n = 1:num_files
        [Gin(:,:,n),tin,fin] = specgram(input(n,:),fs,ns);
        [Gout(:,:,n),t,f] = specgram(output(n,:),fs,ns);
    end
    
else
    Gin = zeros(floor(N/ns), floor(ns/4), length(IDnum));
    Gout = zeros(floor(N/ns), floor(ns/4), length(IDnum));

    for n = 1:length(IDnum)
        [Gin(:,:,n),tin,fin] = specgram(input(n,:),fs,ns);
        [Gout(:,:,n),t,f] = specgram(output(n,:),fs,ns);
    end
end

disp('Autospectra Calculated')

%% calibration

if all_files == true
    Cout = zeros(floor(N/ns), floor(ns/4), num_files); 
    
    for n = 1:num_files
        Cout(:,:,n) = UW_Sensitivity(f, Gout(:,:,n), 1, 4038, 4034);
    end
    
else
    Cout = zeros(floor(N/ns), floor(ns/4), length(IDnum));
    
    for n = 1:length(IDnum)
        Cout(:,:,n) = UW_Sensitivity(f, Gout(:,:,n), 1, 4038, 4034);
    end
end
        
disp('Recieved Signal Calibrated')

%% xcorr

if all_files == false
    if XandA == true
        time_avg = zeros(length(IDnum),N);
        time_avg(1,:) = output(1,:);
        for n = 2:length(IDnum)
            [c,lags] = xcorr(output(1,:), output(n,:));
            [~,indmax] = max(c);
            time_avg(n,:) = circshift(output(n,:),lags(indmax));
        end
        
        [Gavg,tavg,favg] = specgram(mean(time_avg),fs,ns);
        
        disp('Cross-correlated and Time-averaged')
    end
end



%% plotting

if all_files == false
    if XandA == true
        prompt = {'Do you want to see waveforms? Y/N [Y]: ','Do you want to see spectrograms? Y/N [Y]: '...
            'Do you want to see the time averaged waveform? Y/N [Y]: ','Do you want to see the time averaged spectrogram? Y/N [Y]: '};
        dlgtitle = 'Input';
        dims = [1 50];
        User_Input = inputdlg(prompt,dlgtitle,dims)
    else
        prompt = {'Do you want to see waveforms? Y/N [Y]: ','Do you want to see spectrograms? Y/N [Y]: '};
        dlgtitle = 'Input';
        dims = [1 35];
        User_Input = inputdlg(prompt,dlgtitle,dims)
    end
else
    prompt = {'Do you want to see waveforms? Y/N [Y]: ','Do you want to see spectrograms? Y/N [Y]: '};
    dlgtitle = 'Input';
    dims = [1 35];
    User_Input = inputdlg(prompt,dlgtitle,dims)
end

if string(User_Input(1)) == 'Y'
    time = (0:1/N:1-1/N);
    if all_files == true
        
        for n = 1:num_files

            figure()
            plot((time*t_rec).', input(n,:))
            title(sprintf('generated signal ID:%02d',n))
            xlabel('Time (s)')
            ylabel('Amplitude (mV)')
            

            figure()
            plot((time*t_rec).', output(n,:))
            title(sprintf('received signal ID:%02d',n))
            xlabel('Time (s)')
            ylabel('Amplitude (mV)')
        end
    else
      for n = 1:length(IDnum)

            figure()
            plot((time*t_rec).', input(n,:))
            title(sprintf('generated signal ID:%02d',IDnum(n)))
            xlabel('Time (s)')
            ylabel('Amplitude (mV)')
            

            figure()
            plot((time*t_rec).', output(n,:))
            title(sprintf('received signal ID:%02d',IDnum(n)))
            xlabel('Time (s)')
            ylabel('Amplitude (mV)')
      end  
    end
end
            

if string(User_Input(2)) == 'Y'
    if all_files == true
        
        for n = 1:num_files

            figure()
            pcolor(t,f/1000,(20*log10(Gin(:,:,n)/1e-6)).')
            shading interp;
            colorbar;
            colormap('jet')
            title(sprintf('generated signal ID:%02d',n))
            xlabel('time (s)')
            ylabel('frequency(kHz)')

            figure()
            pcolor(t,f/1000,(20*log10(real(Cout(:,:,n))/1e-6)).')
            shading interp;
            colorbar;
            colormap('jet')
            title(sprintf('received signal ID:%02d',n))
            xlabel('time (s)')
            ylabel('frequency(kHz)')
        end
    else
        for n = 1:length(IDnum)
            figure()
            pcolor(t,f/1000,(20*log10(Gin(:,:,n)/1e-6)).')
            shading interp;
            colorbar;
            colormap('jet')
            title(sprintf('generated signal ID:%02d',IDnum(n)))
            xlabel('time (s)')
            ylabel('frequency(kHz)')

            figure()
            pcolor(t,f/1000,(20*log10(real(Cout(:,:,n))/1e-6)).')
            shading interp;
            colorbar;
            colormap('jet')
            title(sprintf('received signal ID:%02d',IDnum(n)))
            xlabel('time (s)')
            ylabel('frequency(kHz)')
        end
    end
end

if all_files == false
    if XandA == true

        if string(User_Input(3)) == 'Y'
            time = (0:1/N:1-1/N);
            figure()
            plot(time*t_rec, mean(time_avg))
            title('time averaged received signal')
            xlabel('Time (s)')
            ylabel('Amplitude (mV)')
        end


        if string(User_Input(4)) == 'Y'
            figure()
            pcolor(tavg,favg/1000,(20*log10(Gavg/1e-6)).')
            shading interp;
            colorbar;
            %caxis([0 100])
            colormap('jet')
            title('time-averaged received signal')
            xlabel('time (s)')
            ylabel('frequency(kHz)')
        end
    end
end
