%% Analog Communications MATLAB Assignment
%% Prepared by:
    %% Ahmed Aly 19015292
    %% Ahmed Sherif 19015255
    %% Yahia Walid 19016891
    %% Youssef Mohamed 19016941
    %% Ali Emad 19016028
%% Introduction

% Clearing Command Window and Shutting all Sounds up
close all;
clc 
clear sound;

%% Prologue: Reading, Filtering, Sketching and Playing

% Introductory Statement
fprintf('Welcome to Project Analogus.\n')
fprintf('\n')

% Reading the Audio File and Storing it in inputAudio
%path = input('Enter the Full Path of the File Including Extension with Quotations ');
%[inputAudio, Fs] =  audioread(path);
[inputAudio, Fs] =  audioread('eric.wav');

% Saving a Single Channel of the Audio File
unfilteredAudio_t = inputAudio(:, 1);

% Converting unfilteredAudio_t from Time to Frequency Domain
unfilteredAudio_f = fftshift(fft(unfilteredAudio_t));

% Defining the Range of Frequencies for unfilteredAudio_f
freqRange_afterResample_dsbsc = linspace(-Fs/2,Fs/2,length(unfilteredAudio_t));

% Defining the Range of Time for unfilteredAudio_t
timeRange_prologue = linspace(0,length(unfilteredAudio_t)/Fs,length(unfilteredAudio_t));
timeRange_prologue = timeRange_prologue';

% Storing both Magnitude and Phase of unfilteredAudio_f Separately
unfilteredAudio_f_mag = abs(unfilteredAudio_f);
unfilteredAudio_f_phase = angle(unfilteredAudio_f)*(180/pi); 

% Plotting unfilteredAudio_t, unfilteredAudio_f_mag, and unfilteredAudio_f_phase
figure
plot(timeRange_prologue,unfilteredAudio_t)
title('Unfiltered Signal - Time Domain')
saveas(gcf,'figures\Exp1\Unfiltered Signal - Time Domain.png')
figure
plot(freqRange_afterResample_dsbsc,unfilteredAudio_f_mag)
title('Unfiltered Signal Magnitude - Frequency Domain')
saveas(gcf,'figures\Exp1\Unfiltered Signal Magnitude - Frequency Domain.png')
figure
plot(freqRange_afterResample_dsbsc,unfilteredAudio_f_phase)
title('Unfiltered Signal Phase - Frequency Domain')
saveas(gcf,'figures\Exp1\Unfiltered Signal Phase - Frequency Domain.png')

% Playing unfilteredAudio_t
fprintf('Playing the Input Audio Signal\n')
sound(unfilteredAudio_t, Fs);
disp("Press Any Key to Stop Sound and Resume the Program")
pause();
clear sound;

% ILPF with Cutoff Frequency 4kHz
filter_freq = 4000;
N_samples = length(unfilteredAudio_t);
n_samples_per_Hz = N_samples/Fs;
right_band = round((Fs/2-filter_freq)*n_samples_per_Hz);
left_band = (N_samples-right_band+1);
filteredAudio_f = unfilteredAudio_f;
filteredAudio_f([1:right_band left_band:N_samples]) = 0;

% Converting the Real Part of filteredAudio_f to Time Domain
filteredAudio_t = real(ifft(ifftshift(filteredAudio_f)));

% Defining the Range of Time and Frequency for filteredAudio_t and unfilteredAudio_f Respectively
timeRange_prologue = linspace(0,length(filteredAudio_t)/Fs,length(filteredAudio_t));
timeRange_prologue = timeRange_prologue';
freqRange_afterResample_dsbsc = linspace(-Fs/2, Fs/2, N_samples);

% Plotting Filtered Audio Signal in both Time and Frequency Domain
figure
plot(timeRange_prologue,filteredAudio_t)
title('Filtered Signal - Time Domain')
saveas(gcf,'figures\Exp1\Filtered Signal - Time Domain.png')
figure
plot(freqRange_afterResample_dsbsc,abs(filteredAudio_f))
title('Filtered Signal - Frequency Domain')
saveas(gcf,'figures\Exp1\Filtered Signal - Frequency Domain.png')

% Playing Filtered Audio Signal unfilteredAudio_t
fprintf('Playing Audio Signal After Filtering Frequencies Higher than 4kHz\n')
sound(filteredAudio_t, Fs);
disp("Press Any Key to Stop Sound and Resume the Program")
pause();
clear sound;
fprintf('\n')
fprintf('\n')

%%%%%%%%%%%End of Prologue%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Transmission Parameters

% Defining Carrier Frequency Fc = 100kHz, m = Vmax/A = 0.5, Vmax = 0.1961 --> A = 0.3922 
fc = 100000;
fs = fc * 5;
Vmax = max(filteredAudio_t);
m = 0.5;
A = Vmax / m;

downsampleFactor = round(fs / Fs);

% Resampling Signal to Match New Sampling Frequency
message = resample(filteredAudio_t,fs,Fs);

% Redefining Time Range After Resampling
timeRange_afterResample = linspace(0,length(message)/fs,length(message));

% Plotting message
figure
plot(timeRange_afterResample,message)
title('Message After Resampling')
saveas(gcf,'figures\Exp1\Message After Resampling.png')

%% Experiment 1: DSB Modulation
fprintf('Experiment 01: DSB Modulation\n')

%%% Transmission
% Defining Carrier ### Transpose (') is ESSENTIAL to Calculate s(t) Correctly
carrier = (cos(2 * pi * fc * timeRange_afterResample)');

%% For DSB-SC
fprintf('DSB-SC\n')

% s(t) = Acos(2*pi*fc*t) * m(t)
st_dsbsc = message .* (A * carrier);

% Plotting the Modulated Signal
figure
plot(timeRange_afterResample,st_dsbsc)
title('DSB-SC: Modulated Signal in Time Domain')
saveas(gcf,'figures\Exp1\DSB-SC - Modulated Signal in Time Domain.png')

% S(f) = FT{s(t)}
Sf_dsbsc = fftshift(fft(st_dsbsc));

% Defining the Range of Frequency for Sf_dsbsc
freqRange_afterResample_dsbsc = linspace(-fs/2,fs/2,length(st_dsbsc));

% Storing both Magnitude and Phase of Sf_dsbsc Separately
Sf_dsbsc_mag = abs(Sf_dsbsc);
Sf_dsbsc_phase = angle(Sf_dsbsc)*(180/pi);

% Plotting Sf_dsbsc
figure
plot(freqRange_afterResample_dsbsc,Sf_dsbsc_mag)
title('DSB-SC: Modulated Signal in Frequency Domain')
saveas(gcf,'figures\Exp1\DSB-SC - Modulated Signal in Frequency Domain.png')

%%% Receiver (Coherent)
% While Loop to Keep Asking the User for Parameters Regarding Coherent Detection - DSB-SC
flag = 0;
while (flag == 0)
    
    % Chosing the Error in Coherent Receiver
    fprintf('Choose the Setting for the Coherent Receiver in DSB-SC\n')
    fprintf('1. No Noise.\n')
    fprintf('2. Account for SNR\n')
    fprintf('3. Account for Frequency Error = +0.1kHz\n')
    fprintf('4. Account for Phase Error = pi/9\n')
    fprintf('5. Proceed to Next Section\n')
    dsbscChoice_coh = input('Enter the Number Corresponding to your Choice ');
    fprintf('\n')
    
    % Operating Based on Choice
    switch (dsbscChoice_coh)
        
        case 1 % no noise
            % Using Coherent Detector -> Multiply the Message by the Carrier
            receivedMessage_dsbsc_coh_noNoise = carrier .* st_dsbsc;

            % Defining the Used 5th Order Butterworth LPF and Implementing it 
            [b, a] = butter(10, 2*fc/fs); % b -> num, a -> denum
            receivedFilteredMessage_dsbsc_coh_noNoise = filter(b,a,receivedMessage_dsbsc_coh_noNoise) * 5; % Multiply by 5 to Compensate Dividing the Amplitude by 5 in Filter

            % Downsampling the Message to Play it (Each 10 Samples in Resampled Signal Translates to 1 Sample in Original Signal)
            receivedAudio_dsbsc_coh_noNoise = downsample(receivedFilteredMessage_dsbsc_coh_noNoise, downsampleFactor);

            % Redefining the Range of Time
            timeRange_dsbsc_coh_noNoise = linspace(0,length(receivedAudio_dsbsc_coh_noNoise)/fs,length(receivedAudio_dsbsc_coh_noNoise));
           
            % Plotting the Received Audio Signal
            figure
            plot(timeRange_dsbsc_coh_noNoise,receivedAudio_dsbsc_coh_noNoise) 
            title('DSB-SC: Demodulated Filtered Signal in Time Domain (Coherent - No Noise)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal in Time Domain (Coherent - No Noise).png')
            
            % Convert to Frequency and plot
            receivedAudio_dsbsc_coh_noNoise_f = fftshift(fft(receivedAudio_dsbsc_coh_noNoise));
            receivedAudio_dsbsc_coh_noNoise_f_mag = abs(receivedAudio_dsbsc_coh_noNoise_f);
            receivedAudio_dsbsc_coh_noNoise_f_phase = angle(receivedAudio_dsbsc_coh_noNoise_f) * (180/pi);
            freqRange_dsbsc_coh = linspace(0,length(receivedAudio_dsbsc_coh_noNoise)/fs,length(receivedAudio_dsbsc_coh_noNoise));
            figure
            plot(freqRange_dsbsc_coh,receivedAudio_dsbsc_coh_noNoise_f_mag) 
            title('DSB-SC: Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent - No Noise)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent - No Noise).png')
            figure
            plot(freqRange_dsbsc_coh,receivedAudio_dsbsc_coh_noNoise_f_phase) 
            title('DSB-SC: Demodulated Filtered Signal Phase in Frequency Domain (Coherent - No Noise)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Phase in Frequency Domain (Coherent - No Noise).png')


            % Playing The Received Audio Signal
            fprintf('Playing Received Audio Signal - DSB-SC (Coherent - No Noise)\n')
            sound(receivedAudio_dsbsc_coh_noNoise, Fs);
            disp("Press Any Key to Stop Sound and Resume the Program")
            pause();
            clear sound;
            
        case 2 % SNR
            % SNR Value
            SNR = input('Enter the SNR Value ');
            
            % Applying SNR on st_dsbsc
            st_dsbsc_snr = awgn(st_dsbsc,SNR);
            
            % Using Coherent Detector -> Multiply the Message by the Carrier
            receivedMessage_dsbsc_coh_SNR = carrier .* st_dsbsc_snr;

            % Defining the Used 5th Order Butterworth LPF and Implementing it 
            [b, a] = butter(10, 2*fc/fs); % b -> num, a -> denum
            receivedFilteredMessage_dsbsc_coh_SNR = filter(b,a,receivedMessage_dsbsc_coh_SNR) * 5; % Multiply by 5 to Compensate Dividing the Amplitude by 5 in Filter

            % Downsampling the Message to Play it (Each 10 Samples in Resampled Signal Translates to 1 Sample in Original Signal)
            receivedAudio_dsbsc_coh_SNR = downsample(receivedFilteredMessage_dsbsc_coh_SNR, downsampleFactor);

            % Redefining the Range of Time
            timeRange_dsbsc_coh_SNR = linspace(0,length(receivedAudio_dsbsc_coh_SNR)/fs,length(receivedAudio_dsbsc_coh_SNR));

            % Plotting the Received Audio Signal
            figure
            plot(timeRange_dsbsc_coh_SNR,receivedAudio_dsbsc_coh_SNR) 
            title('DSB-SC: Demodulated Filtered Signal in Time Domain (Coherent - With SNR)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal in Time Domain (Coherent - With SNR).png')

            
            % Convert to Frequency and plot
            receivedAudio_dsbsc_coh_SNR_f = fftshift(fft(receivedAudio_dsbsc_coh_SNR));
            receivedAudio_dsbsc_coh_SNR_f_mag = abs(receivedAudio_dsbsc_coh_SNR_f);
            receivedAudio_dsbsc_coh_SNR_f_phase = angle(receivedAudio_dsbsc_coh_SNR_f) * (180/pi);
            freqRange_dsbsc_coh = linspace(0,length(receivedAudio_dsbsc_coh_SNR)/fs,length(receivedAudio_dsbsc_coh_SNR));
            figure
            plot(freqRange_dsbsc_coh,receivedAudio_dsbsc_coh_SNR_f_mag) 
            title('DSB-SC: Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent - With SNR)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent - With SNR).png')

            figure
            plot(freqRange_dsbsc_coh,receivedAudio_dsbsc_coh_SNR_f_phase) 
            title('DSB-SC: Demodulated Filtered Signal Phase in Frequency Domain (Coherent - With SNR)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Phase in Frequency Domain (Coherent - With SNR).png')


            % Playing The Received Audio Signal
            fprintf('Playing Received Audio Signal - DSB-SC (Coherent - With SNR)\n')
            sound(receivedAudio_dsbsc_coh_SNR, Fs);
            disp("Press Any Key to Stop Sound and Resume the Program")
            pause();
            clear sound;
        
        case 3 % Frequency Error
            % Frequency Error
            freqError = 100;
            
            % Carrier with Frequency Error
            carrier_freqError = (cos(2 * pi * (fc + freqError) * timeRange_afterResample)');
            
            % Using Coherent Detector -> Multiply the Message by the Carrier
            receivedMessage_dsbsc_coh_freqError = carrier_freqError .* st_dsbsc;

            % Defining the Used 5th Order Butterworth LPF and Implementing it 
            [b, a] = butter(10, 2*fc/fs); % b -> num, a -> denum
            receivedFilteredMessage_dsbsc_coh_freqError = filter(b,a,receivedMessage_dsbsc_coh_freqError); % Multiply by 5 to Compensate Dividing the Amplitude by 5 in Filter

            % Downsampling the Message to Play it (Each 10 Samples in Resampled Signal Translates to 1 Sample in Original Signal)
            receivedAudio_dsbsc_coh_freqError = downsample(receivedFilteredMessage_dsbsc_coh_freqError, downsampleFactor);

            % Redefining the Range of Time
            timeRange_dsbsc_coh_freqError = linspace(0,length(receivedAudio_dsbsc_coh_freqError)/fs,length(receivedAudio_dsbsc_coh_freqError));

            % Plotting the Received Audio Signal
            figure
            plot(timeRange_dsbsc_coh_freqError,receivedAudio_dsbsc_coh_freqError) 
            title('DSB-SC: Demodulated Filtered Signal in Time Domain (Coherent - With Frequency Error = 0.1kHz)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal in Time Domain (Coherent - With Frequency Error = 0.1kHz).png')

            
            % Convert to Frequency and plot
            receivedAudio_dsbsc_coh_freqError_f = fftshift(fft(receivedAudio_dsbsc_coh_freqError));
            receivedAudio_dsbsc_coh_freqError_f_mag = abs(receivedAudio_dsbsc_coh_freqError_f);
            receivedAudio_dsbsc_coh_freqError_f_phase = angle(receivedAudio_dsbsc_coh_freqError_f) * (180/pi);
            freqRange_dsbsc_coh = linspace(0,length(receivedAudio_dsbsc_coh_freqError)/fs,length(receivedAudio_dsbsc_coh_freqError));
            figure
            plot(freqRange_dsbsc_coh,receivedAudio_dsbsc_coh_freqError_f_mag) 
            title('DSB-SC: Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent - With Frequency Error = 0.1kHz)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent - With Frequency Error = 0.1kHz).png')

            figure
            plot(freqRange_dsbsc_coh,receivedAudio_dsbsc_coh_freqError_f_phase) 
            title('DSB-SC: Demodulated Filtered Signal Phase in Frequency Domain (Coherent - With Frequency Error = 0.1kHz)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Phase in Frequency Domain (Coherent - With Frequency Error = 0.1kHz).png')


            % Playing The Received Audio Signal
            fprintf('Playing Received Audio Signal - DSB-SC (Coherent - With Frequency Error)\n')
            sound(receivedAudio_dsbsc_coh_freqError, Fs);
            disp("Press Any Key to Stop Sound and Resume the Program")
            pause();
            clear sound;
        
        case 4 % Phase Error
            % Phase Error
            phaseError = 20 * pi / 180;
            
            % Carrier with Phase Error
            carrier_phaseError = (cos((2 * pi * fc * timeRange_afterResample) + phaseError)');
            
            % Using Coherent Detector -> Multiply the Message by the Carrier
            receivedMessage_dsbsc_coh_phaseError = carrier_phaseError .* st_dsbsc;

            % Defining the Used 5th Order Butterworth LPF and Implementing it 
            [b, a] = butter(10, 2*fc/fs); % b -> num, a -> denum
            receivedFilteredMessage_dsbsc_coh_phaseError = filter(b,a,receivedMessage_dsbsc_coh_phaseError) * 5; % Multiply by 5 to Compensate Dividing the Amplitude by 5 in Filter

            % Downsampling the Message to Play it (Each 10 Samples in Resampled Signal Translates to 1 Sample in Original Signal)
            receivedAudio_dsbsc_coh_phaseError = downsample(receivedFilteredMessage_dsbsc_coh_phaseError, downsampleFactor);

            % Redefining the Range of Time
            timeRange_dsbsc_coh_phaseError = linspace(0,length(receivedAudio_dsbsc_coh_phaseError)/fs,length(receivedAudio_dsbsc_coh_phaseError));

            % Plotting the Received Audio Signal
            figure
            plot(timeRange_dsbsc_coh_phaseError,receivedAudio_dsbsc_coh_phaseError) 
            title('DSB-SC: Demodulated Filtered Signal in Time Domain (Coherent - With Phase Error = pi/9)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal in Time Domain (Coherent - With Phase Error = pi%9).png')

            
            % Convert to Frequency and plot
            receivedAudio_dsbsc_coh_phaseError_f = fftshift(fft(receivedAudio_dsbsc_coh_phaseError));
            receivedAudio_dsbsc_coh_phaseError_f_mag = abs(receivedAudio_dsbsc_coh_phaseError_f);
            receivedAudio_dsbsc_coh_phaseError_f_phase = angle(receivedAudio_dsbsc_coh_phaseError_f) * (180/pi);
            freqRange_dsbsc_coh = linspace(0,length(receivedAudio_dsbsc_coh_phaseError)/fs,length(receivedAudio_dsbsc_coh_phaseError));
            figure
            plot(freqRange_dsbsc_coh,receivedAudio_dsbsc_coh_phaseError_f_mag) 
            title('DSB-SC: Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent - With Phase Error = pi/9)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent - With Phase Error = pi%9).png')

            figure
            plot(freqRange_dsbsc_coh,receivedAudio_dsbsc_coh_phaseError_f_phase) 
            title('DSB-SC: Demodulated Filtered Signal Phase in Frequency Domain (Coherent - With Phase Error = pi/9)')
            saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Phase in Frequency Domain (Coherent - With Phase Error = pi%9).png')


            % Playing The Received Audio Signal
            fprintf('Playing Received Audio Signal - DSB-SC (Coherent - With Phase Error)\n')
            sound(receivedAudio_dsbsc_coh_phaseError, Fs);
            disp("Press Any Key to Stop Sound and Resume the Program")
            pause();
            clear sound;
            
        case 5 % Continue Program
            flag = 1;
   
        otherwise % Idiotproofing
            fprintf('Invalid Input! Try Again\n')
    end
    
end

%%% Receiver (Envelope)
% Using Envelope Detector
envelopeReceiver_dsbsc = abs(hilbert(st_dsbsc));

% Downsampling the Message to Play it
receivedAudio_dsbsc_env = downsample(envelopeReceiver_dsbsc, downsampleFactor);

% Redefining the Range of Time
timeRange_dsbsc_env = linspace(0,length(receivedAudio_dsbsc_env)/fs,length(receivedAudio_dsbsc_env));

% Plotting the Received Message (Distorted)
figure
plot(timeRange_dsbsc_env,receivedAudio_dsbsc_env) 
title('DSB-SC: Demodulated Filtered Signal in Time Domain (Envelope)')
saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal in Time Domain (Envelope).png')


% Convert to Frequency and Plot
receivedAudio_dsbsc_env_f = fftshift(fft(receivedAudio_dsbsc_env));
receivedAudio_dsbsc_env_f_mag = abs(receivedAudio_dsbsc_env_f);
receivedAudio_dsbsc_env_f_phase = angle(receivedAudio_dsbsc_env_f) * (180/pi);
freqRange_dsbsc_env = linspace(0,length(receivedAudio_dsbsc_env)/fs,length(receivedAudio_dsbsc_env));
figure
plot(freqRange_dsbsc_env,receivedAudio_dsbsc_env_f_mag) 
title('DSB-SC: Demodulated Filtered Signal Magnitude in Frequency Domain (Envelope)')
saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Magnitude in Frequency Domain (Envelope).png')

figure
plot(freqRange_dsbsc_env,receivedAudio_dsbsc_env_f_phase) 
title('DSB-SC: Demodulated Filtered Signal Phase in Frequency Domain (Envelope)')
saveas(gcf,'figures\Exp1\DSB-SC - Demodulated Filtered Signal Phase in Frequency Domain (Envelope).png')

% Not gonna play it cuz it's distorted
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For DSB-TC
fprintf('DSB-TC\n')

% s(t)
st_dsbtc = ((A .* carrier) + (message .* carrier));

% Plotting the Modulated Signal
figure
plot(timeRange_afterResample,st_dsbtc)
title('DSB-TC: Modulated Signal in Time Domain')
saveas(gcf,'figures\Exp1\DSB-TC - Modulated Signal in Time Domain.png')

% S(f)
Sf_dsbtc = fftshift(fft(st_dsbtc));

% Defining the Range of Frequencies for Sf_tc
freqRange_afterResample_dsbtc = linspace(-fs/2,fs/2,length(st_dsbtc));

% Storing both Magnitude and Phase of Sf_dsbtc Separately
Sf_dsbtc_mag = abs(Sf_dsbtc);
Sf_dsbtc_phase = angle(Sf_dsbtc)*(180/pi);

% Plotting Sf_dsbtc_mag
figure
plot(freqRange_afterResample_dsbtc,Sf_dsbtc_mag)
title('DSB-TC: Modulated Signal in Frequency Domain')
saveas(gcf,'figures\Exp1\DSB-TC - Modulated Signal in Frequency Domain.png')


%%% Receiving (Coherent)
% Using Coherent Detector -> Multiply the Message by the Carrier
receivedMessage_dsbtc_coh = carrier .* st_dsbtc;

% Defining the Used 5th Order Butterworth LPF and Implementing it
[b, a] = butter(5, 2*fc/fs); % b -> num, a -> denum
receivedFilteredMessage_dsbtc_coh = filter(b,a,receivedMessage_dsbtc_coh); 

% Multiply by 2 to Compensate Dividing the Amplitude by 2 in Filter
receivedFilteredMessage_dsbtc_coh = receivedFilteredMessage_dsbtc_coh * 2;

% DC Blocking
receivedFilteredMessage_dsbtc_coh = receivedFilteredMessage_dsbtc_coh - A;

% Downsampling the Message to Play it
receivedAudio_dsbtc_coh = downsample(receivedFilteredMessage_dsbtc_coh, 10);

% Redefining the Range of Time
timeRange_dsbtc_coh = linspace(0,length(receivedAudio_dsbtc_coh)/fs,length(receivedAudio_dsbtc_coh));

% Plotting the Received Signal
figure
plot(timeRange_dsbtc_coh,receivedAudio_dsbtc_coh) % multiply by 2 to compensate the for dividing the amplitude by 2
title('DSB-TC: Demodulated Filtered Signal in Time Domain (Coherent)')
saveas(gcf,'figures\Exp1\DSB-TC - Demodulated Filtered Signal in Time Domain (Coherent).png')


% Convert to Frequency and Plot
receivedAudio_dsbtc_coh_f = fftshift(fft(receivedAudio_dsbtc_coh));
receivedAudio_dsbtc_coh_f_mag = abs(receivedAudio_dsbtc_coh_f);
receivedAudio_dsbtc_coh_f_phase = angle(receivedAudio_dsbtc_coh_f) * (180/pi);
freqRange_dsbtc_coh = linspace(0,length(receivedAudio_dsbtc_coh)/fs,length(receivedAudio_dsbtc_coh));
figure
plot(freqRange_dsbtc_coh,receivedAudio_dsbtc_coh_f_mag) 
title('DSB-TC: Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent)')
saveas(gcf,'figures\Exp1\DSB-TC - Demodulated Filtered Signal Magnitude in Frequency Domain (Coherent).png')

figure
plot(freqRange_dsbtc_coh,receivedAudio_dsbtc_coh_f_phase) 
title('DSB-TC: Demodulated Filtered Signal Phase in Frequency Domain (Coherent)')
saveas(gcf,'figures\Exp1\DSB-TC - Demodulated Filtered Signal Phase in Frequency Domain (Coherent).png')

% Playing The Received Audio Signal
fprintf('Playing Received Audio Signal - DSB-TC (Coherent)\n')
sound(receivedAudio_dsbtc_coh, Fs);
disp("Press Any Key to Stop Sound and Resume the Program")
pause();
clear sound;

%%%%%%%%%%%%%

% Receiver (Envelope)
% Using Envelope Detector
envelopeReceiver_dsbtc = abs(hilbert(st_dsbtc));

% Downsampling the Message to Play it
receivedAudio_dsbtc_env = downsample(envelopeReceiver_dsbtc, 10);

% DC Blocking
receivedAudio_dsbtc_env = receivedAudio_dsbtc_env - A;

% Redefining the Range of Time
timeRange_dsbtc_env = linspace(0,length(receivedAudio_dsbtc_env)/fs,length(receivedAudio_dsbtc_env));

% Plotting the Received Message (Good)
figure
plot(timeRange_dsbtc_env,receivedAudio_dsbtc_env) 
title('DSB-TC: Demodulated Filtered Signal in Time Domain (Envelope)')
saveas(gcf,'figures\Exp1\DSB-TC - Demodulated Filtered Signal in Time Domain (Envelope).png')

% Convert to Frequency and Plot
receivedAudio_dsbtc_env_f = fftshift(fft(receivedAudio_dsbtc_env));
receivedAudio_dsbtc_env_f_mag = abs(receivedAudio_dsbtc_env_f);
receivedAudio_dsbtc_env_f_phase = angle(receivedAudio_dsbtc_env_f) * (180/pi);
freqRange_dsbtc_env = linspace(0,length(receivedAudio_dsbtc_env)/fs,length(receivedAudio_dsbtc_env));
figure
plot(freqRange_dsbtc_env,receivedAudio_dsbtc_env_f_mag)
title('DSB-TC: Demodulated Filtered Signal Magnitude in Frequency Domain (Envelope)')
saveas(gcf,'figures\Exp1\DSB-TC - Demodulated Filtered Signal Magnitude in Frequency Domain (Envelope).png')

figure
plot(freqRange_dsbtc_env,receivedAudio_dsbtc_env_f_phase) 
title('DSB-TC: Demodulated Filtered Signal Phase in Frequency Domain (Envelope)')
saveas(gcf,'figures\Exp1\DSB-TC - Demodulated Filtered Signal Phase in Frequency Domain (Envelope).png')

% Playing The Received Audio Signal
fprintf('Playing Received Audio Signal - DSB-TC (Envelope)\n')
sound(receivedAudio_dsbtc_env, Fs);
disp("Press Any Key to Stop Sound and Resume the Program")
pause();
clear sound;
fprintf('\n')
fprintf('\n')


%%%%%%%%%%%End of DSB Modulation%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Experiment 02: SSB Modulation
fprintf('Experiment 02: SSB Modulation\n')

%%% Transmission
% Defining Carrier ### Transpose (') is ESSENTIAL to Calculate s(t) Correctly
carrier = (cos(2 * pi * fc * timeRange_afterResample)');
carrier_hat = hilbert(carrier);
message_hat = hilbert(message);

%% SSB-SC
fprintf('SSB-SC\n')

% Add to Keep LSB -> s(t) = m(t).carrier + m~(t).carrier_sin
st_ssbsc = real((message .* carrier) + (message_hat .* carrier_hat));

% Plotting the Modulated Signal
figure
plot(timeRange_afterResample,st_ssbsc)
title('SSB-SC: Modulated Signal in Time Domain')

% S(f) = FT{s(t)}
receivedMessage_ssbsc_coh_SNR_F = fftshift(fft(st_ssbsc));

% Defining the Range of Frequencies for Sf_ssbsc
freqRange_afterResample_ssbsc = linspace(-fs/2,fs/2,length(st_ssbsc));

% Storing both Magnitude and Phase of Sf_ssbsc Separately
Sf_ssbsc_mag = abs(receivedMessage_ssbsc_coh_SNR_F);
Sf_ssbsc_phase = angle(receivedMessage_ssbsc_coh_SNR_F)*(180/pi);

% Plotting Sf_ssbsc
figure
plot(freqRange_afterResample_ssbsc,Sf_ssbsc_mag)
title('SSB-SC: Modulated Signal in Frequency Domain')

%%% Receiving (Coherent)
% While Loop to Keep Asking the User for Parameters Regarding Coherent Detection - DSB-SC
flag = 0;
while (flag == 0)
    % Chosing the Error in Coherent Receiver
    fprintf('Choose the Filter for the Coherent Receiver in SSB-SC\n')
    fprintf('1. Butterworth Filter\n')
    fprintf('2. Ideal Low Pass Filter\n')
    fprintf('3. Proceed to Next Section\n')
    ssbscChoice_coh = input('Enter the Number Corresponding to your Choice ');
    fprintf('\n')
    
    switch (ssbscChoice_coh) 
        
        case 1 % Butterworth
            % Multiply the Message by the Carrier
            receivedMessage_ssbsc_coh_butter = carrier .* st_ssbsc;

            % Defining the Used 4th Order Butterworth LPF and Implementing it
            [b, a] = butter(4, 2*fc/fs); % b -> num, a -> denum

            % Multiply by 2 to Compensate for Dividing the Amplitude by 2
            % Applying the Filter
            receivedFilteredMessage_ssbsc_coh_butter = filter(b,a,receivedMessage_ssbsc_coh_butter);
               
            % Downsampling the Message to Play it
            receivedFilteredMessage_ssbsc_coh_butter = downsample(receivedFilteredMessage_ssbsc_coh_butter, downsampleFactor);

            % Redefining the Range of Time
            timeRange_ssbsc_coh_butter = linspace(0,length(receivedFilteredMessage_ssbsc_coh_butter)/fs,length(receivedFilteredMessage_ssbsc_coh_butter));
            
            % Plotting the Received Audio Signal
            figure
            plot(timeRange_ssbsc_coh_butter,receivedFilteredMessage_ssbsc_coh_butter) % multiply by 5 to compensate the for dividing the amplitude by 5
            title('SSB-SC: Demodulated Filtered Signal in Time Domain (Butterworth)')
            
            % Convert to Frequency and Plot
            receivedFilteredMessage_ssbsc_coh_butter_f = fftshift(fft(receivedFilteredMessage_ssbsc_coh_butter));
            receivedFilteredMessage_ssbsc_coh_butter_f_mag = abs(receivedFilteredMessage_ssbsc_coh_butter_f);
            receivedFilteredMessage_ssbsc_coh_butter_f_phase = angle(receivedFilteredMessage_ssbsc_coh_butter_f) * (180/pi);
            freqRange_ssbsc_coh = linspace(0,length(receivedFilteredMessage_ssbsc_coh_butter)/fs,length(receivedFilteredMessage_ssbsc_coh_butter));
            figure
            plot(freqRange_ssbsc_coh,receivedFilteredMessage_ssbsc_coh_butter_f_mag)
            title('SSB-SC: Demodulated Filtered Signal Magnitude in Frequency Domain (Butterworth)')
            figure
            plot(freqRange_ssbsc_coh,receivedFilteredMessage_ssbsc_coh_butter_f_phase) 
            title('SSB-SC: Demodulated Filtered Signal Phase in Frequency Domain (Butterworth)')

            % Playing The Received Audio Signal
            fprintf('Playing Received Audio Signal - SSB-SC (Butterworth)\n')
            sound(receivedFilteredMessage_ssbsc_coh_butter, Fs);
            disp("Press Any Key to Stop Sound and Resume the Program")
            pause();
            clear sound;
            
            
        case 2 %ILPF
            % ILPF with Cutoff Frequency 100kHz + ADD SNR to received signal ##
            filter_freq = 100000;
            N = length(receivedMessage_ssbsc_coh_SNR_F);
            n = N/Fs;
            right_band = round((Fs/2-filter_freq)*n);
            left_band = (N-right_band+1);
            
            % SNR Value
            SNR = input('Enter the SNR Value ');
            
            % Applying SNR on st_dsbsc
            st_ssbsc_snr = awgn(st_ssbsc,SNR);
            
            % Multiply the Message by the Carrier
            receivedMessage_ssbsc_coh_SNR = carrier .* st_ssbsc_snr;
            
            % Implementing the Used ILPF 
            receivedMessage_ssbsc_coh_SNR_F = fftshift(fft(receivedMessage_ssbsc_coh_SNR));
            receivedMessage_ssbsc_coh_SNR_F([1:right_band left_band:N]) = 0;
            receivedAudio_ssbsc_coh_SNR = real(ifft(ifftshift(receivedMessage_ssbsc_coh_SNR_F)));

            % Downsampling the Message to Play it
            receivedAudio_ssbsc_coh_SNR = downsample(receivedAudio_ssbsc_coh_SNR, 10);

            % Redefining the Range of Time
            timeRange_ssbsc_coh_SNR = linspace(0,length(receivedAudio_ssbsc_coh_SNR)/fs,length(receivedAudio_ssbsc_coh_SNR));
            
            % Plotting the Received Message
            figure
            plot(timeRange_ssbsc_coh_SNR,receivedAudio_ssbsc_coh_SNR)
            title('SSB-SC: Demodulated Filtered Signal in Time Domain (ILPF)')
            
            % Convert to Frequency and Plot
            receivedAudio_ssbsc_coh_SNR_f = fftshift(fft(receivedAudio_ssbsc_coh_SNR));
            receivedAudio_ssbsc_coh_SNR_f_mag = abs(receivedAudio_ssbsc_coh_SNR_f);
            receivedAudio_ssbsc_coh_SNR_f_phase = angle(receivedAudio_ssbsc_coh_SNR_f) * (180/pi);
            freqRange_ssbsc_coh = linspace(0,length(receivedAudio_ssbsc_coh_SNR)/fs,length(receivedAudio_ssbsc_coh_SNR));
            figure
            plot(freqRange_ssbsc_coh,receivedAudio_ssbsc_coh_SNR_f_mag)
            title('SSB-SC: Demodulated Filtered Signal Magnitude in Frequency Domain (ILPF)')
            figure
            plot(freqRange_ssbsc_coh,receivedAudio_ssbsc_coh_SNR_f_phase) 
            title('SSB-SC: Demodulated Filtered Signal Phase in Frequency Domain (ILPF)')
            
            % Playing The Received Audio Signal
            fprintf('Playing Received Audio Signal - SSB-SC (Coherent - ILPF - With SNR)\n')
            sound(receivedAudio_ssbsc_coh_SNR, Fs);
            disp("Press Any Key to Stop Sound and Resume the Program")
            pause();
            clear sound;
            
        case 3
            flag = 1;
            
        otherwise
            fprintf('Invalid Input! Try Again\n')
    end
end
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% For SSB-TC
fprintf('SSB-TC\n')

% Add to Keep LSB -> s(t) = m(t).carrier + m~(t).carrier_sin
st_ssbtc = (A .* carrier) + (message .* carrier) + (message_hat .* carrier_hat);

% Plotting the Modulated Signal
figure
plot(timeRange_afterResample,st_ssbtc)
title('SSB-TC: Modulated Signal in Time Domain')

% S(f) = FT{s(t)}
Sf_ssbtc = fftshift(fft(st_ssbtc));

% Defining the Range of Frequencies for Sf_ssbtc
freqRange_afterResample_ssbtc = linspace(-fs/2,fs/2,length(st_ssbtc));

% Storing both Magnitude and Phase of Sf_ssbtc Separately
Sf_ssbtc_mag = abs(Sf_ssbtc);
Sf_ssbtc_phase = angle(Sf_ssbtc)*(180/pi);

% Plotting Sf_dsbtc_mag
figure
plot(freqRange_afterResample_ssbtc,Sf_ssbtc_mag)
title('SSB-TC: Modulated Signal in Frequency Domain')


%%% Receiver (Envelope)
% Using Envelope Detector
envelopeReceiver_ssbtc = abs(hilbert(st_ssbtc));

% DC Blocking
envelopeReceiver_ssbtc = envelopeReceiver_ssbtc - A;

% Downsampling the Message to Play it
receivedAudio_ssbtc_env = downsample(envelopeReceiver_ssbtc, 10);

% Redefining the Range of Time
timeRange_ssbtc_env = linspace(0,length(receivedAudio_ssbtc_env)/fs,length(receivedAudio_ssbtc_env));

% Plotting the Received Message
figure
plot(timeRange_ssbtc_env,receivedAudio_ssbtc_env) 
title('SSB-TC: Demodulated Filtered Signal in Time Domain (Envelope)')

% Convert to Frequency and Plot
receivedAudio_ssbtc_env_f = fftshift(fft(receivedAudio_ssbtc_env));
receivedAudio_ssbtc_env_f_mag = abs(receivedAudio_ssbtc_env_f);
receivedAudio_ssbtc_env_f_phase = angle(receivedAudio_ssbtc_env_f) * (180/pi);
freqRange_ssbtc_env = linspace(0,length(receivedAudio_ssbtc_env)/fs,length(receivedAudio_ssbtc_env));
figure
plot(freqRange_ssbtc_env,receivedAudio_ssbtc_env_f_mag)
title('SSB-SC: Demodulated Filtered Signal Magnitude in Frequency Domain (Butterworth)')
figure
plot(freqRange_ssbtc_env,receivedAudio_ssbtc_env_f_phase) 
title('SSB-SC: Demodulated Filtered Signal Phase in Frequency Domain (Butterworth)')

% Playing The Received Audio Signal
fprintf('Playing Received Audio Signal - SSB-TC (Envelope)\n')
sound(receivedAudio_ssbtc_env, Fs);
disp("Press Any Key to Stop Sound and Resume the Program")
pause();
clear sound;
fprintf('\n')
fprintf('\n')


%% Experiment 03: NBFM
fprintf('Experiment 03: NBFM\n')
fprintf('\n')


%%% Transmission
Kf = input('Enter a value for Kf ');
A=input('Enter Gain Value ');
% Defining Carrier ### Transpose (') is ESSENTIAL to Calculate s(t) Correctly
carrier = A.*(cos(2 * pi * fc * timeRange_afterResample)');
carrier_sin = A.*(sin(2 * pi * fc * timeRange_afterResample)');

% Message Integration
message_integration = cumsum(message);
scale=Kf*max(message);
% Calculating Narrow Band Frequency-Modulated Signal
NBFM = ( carrier - ( Kf .*message_integration .* carrier_sin) );

figure
plot(timeRange_afterResample,real(NBFM))
title('NBFM: Modulated Signal in Time Domain')

% Converting to Frequency Domain
NBFM_f = fftshift(fft(real(NBFM)));

% Defining the Range of Frequency for F_NBFM
freqRange_afterResample_NBFM = linspace(-fs/2,fs/2,length(NBFM));

% Storing both Magnitude and Phase of F_NBFM Separately
NBFM_f_mag = abs(NBFM_f);
NBFM_f_phase = angle(NBFM_f) * (180/pi);

% Plotting Sf_dsbsc
figure
plot(freqRange_afterResample_NBFM,NBFM_f_mag)
title('NBFM: Modulated Signal Magnitude in Frequency Domain')

%figure
%plot(freqRange_afterResample_NBFM,NBFM_f_phase)
%title('NBFM: Modulated Signal Phase in Frequency Domain')
%saveas(gcf,fullfile(fname, 'NBFM Modulated Signal Phase in Frequency Domain.png'))

%%% Receiver (Discriminator/Slope Detector)
% Defining Envelope
envelope_NBFM = abs(hilbert(NBFM));
%B = A * message_integration;
%envelope_NBFM = sqrt(A.^2 + B.^2);



receivedMessage_NBFM = zeros(length(NBFM),1);
receivedMessage_NBFM(2:end) = diff(envelope_NBFM);

% Downsampling the Message to Play it
receivedAudio_NBFM = downsample(receivedMessage_NBFM, downsampleFactor);

% Redefining the Range of Time
timeRange_NBFM = linspace(0,length(receivedAudio_NBFM)*10/fs,length(receivedAudio_NBFM));
receivedAudio_NBFM=receivedAudio_NBFM/scale;


%lowpass Butterworth filter

[b,a] = butter(4,2000/(fs/20));
receivedAudio_NBFM = filter(b,a,receivedAudio_NBFM);


% Plotting the Received Message
figure
plot(timeRange_NBFM,receivedAudio_NBFM) 
title('NBFM: Demodulated Received Signal in Time Domain')

% Convert to Frequency and Plot
receivedAudio_NBFM_f = fftshift(fft(receivedAudio_NBFM));
receivedAudio_NBFM_f_mag = abs(receivedAudio_NBFM_f);
receivedAudio_NBFM_f_phase = angle(receivedAudio_NBFM_f);
freqRange_NBFM = linspace(-fs/2,fs/2,length(receivedAudio_NBFM));
figure
plot(freqRange_NBFM,receivedAudio_NBFM_f_mag) 
title('NBFM: Demodulated Received Signal Magnitude in Frequency Domain')

%figure
%plot(freqRange_NBFM,receivedAudio_NBFM_f_phase) 
%title('NBFM: Demodulated Received Signal Phase in Frequency Domain')
%saveas(gcf,fullfile(fname, 'NBFM Demodulated Received Signal Phase in Frequency Domain.png'))

% Playing The Received Audio Signal
fprintf('Playing Received Audio Signal - NBFM \n')
sound(receivedAudio_NBFM, Fs);
disp("Press Any Key to Stop Sound and Resume the Program")
pause();
clear sound;
fprintf('\n')
fprintf('\n')



%%%%%%%%%%%%%%%%End of NBFM%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Epilogue

% Goodbye Message
disp("Thank You for Using Project Analogus")
fprintf('\n')
disp("Adios")
