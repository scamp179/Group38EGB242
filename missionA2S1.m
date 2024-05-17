%% EGB242 Assignment 2, Section 1 %%
% This file is a template for your MATLAB solution to Section 1.
%
% Before starting to write code, generate your data with the ??? as
% described in the assignment task.

%% Initialise workspace
clear all; close all;
load DataA2 audioMultiplexNoisy fs sid;
%% 1.1 
% Compute the time vector using the helper function
t = timevec(0, length(audioMultiplexNoisy)/fs, length(audioMultiplexNoisy));

% Compute the frequency vector using the helper function
f = freqvec(fs, length(audioMultiplexNoisy));

% Compute the FFT of the audio signal
audio_noisy_fft = fftshift(fft(audioMultiplexNoisy)) / fs;

figure;
subplot(2, 1, 1);
plot(t, audioMultiplexNoisy);
title('Time Domain Representation of audioMultiplexNoisy');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(f, abs(audio_noisy_fft));
title('Frequency Domain Representation of audioMultiplexNoisy');
xlabel('Frequency (Hz)');
ylabel('Magnitude');


%% 1.2 modulate and demodulate each audio freq
freq_1 = 8310;
freq_2 = 24090;
freq_3 = 40140;
freq_4 = 56300;
freq_5 = 72110;
% Carrier frequencies identified from the plot
carrierFreqs = [freq_1,freq_2,freq_3,freq_4,freq_5]; % Frequencies in Hz

for idx = 1:length(carrierFreqs)
    % Current carrier frequency
    carrier_freq = carrierFreqs(idx);

    % cutoff frequency
    fc = 2800; % try 3500 to compare
    
    % Demodulate the signal by multiplying it with a cosine wave of the carrier frequency.
    demodulatedSignal = audioMultiplexNoisy .* cos(2*pi*carrier_freq*t);
    
    % Apply a low-pass filter to the demodulated signal to remove high-frequency components
    audioReceived = lowpass(demodulatedSignal, fc, fs);
    
    % Store the demodulated, filtered audio for later use
    audioReceivedCell{idx} = audioReceived;

    audioReceived_shift = fftshift(fft(audioReceived))/fs;
    
    % Plot the filtered (received) audio signal.
    figure;

    % Subplot for the time-domain representation
    subplot(2, 1, 1); % Two rows, one column, first subplot
    plot(t, audioReceived);
    title(sprintf('Time Domain - Carrier %d Hz', carrier_freq));
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(2,1,2);
    plot(f, abs(audioReceived_shift));
    title('Frequency Domain Representation of Received Audio', carrier_freq);
    xlabel('Frequency (HZ)');
    ylabel('Magnitude');


end

%%
%Optional: Listen to each audioReceived after plotting
% for idx = 1:length(audioReceivedCell)/2
%     sound(audioReceivedCell{idx}, fs);
%     pause(length(audioReceivedCell{idx})/fs + 1);  % Play each sound with a pause
% end

%% 1.3 Model the frequency-dependent distortion

% Create an impulse signal
impulse = zeros(1, length(audio_noisy_fft));  % Vector of zeros

impulse(1) = 1;            % Set the first sample to 1

% Transmit the impulse through the channel
output = channel(sid, impulse, fs);

% Frequency response of the impulse response
output_fft = fft(output);  



% Plotting the frequency responses on the same axes
figure;
hold on;
plot(f, real(output_fft), 'r');  % Channel impulse response in blue
plot(f, real(audio_noisy_fft), 'b'); % Multiplexed signal in red
title('Comparison of Frequency Responses');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
legend('Channel Impulse Response', 'Multiplexed Audio Signal');
hold off;
%% 1.4
% we need to find X(F) of noisy audio
noisy_fft = fftshift(fft(audioMultiplexNoisy)); % this in frequency domain Y(F)

%Y(F)/H(F) = X(F)
Denoised = noisy_fft ./output_fft; % = X(F)

Denoised_ifft = ifft(ifftshift(Denoised))* fs; % x(t)



for idx = 1:length(carrierFreqs)
    % Current carrier frequency
    carrier_freq = carrierFreqs(idx);

    % cutoff frequency
    fc = 2000; % try 3500 to compare
    
    % Demodulate the signal by multiplying it with a cosine wave of the carrier frequency.
    demodulated_clean = Denoised_ifft .* cos(2*pi*carrier_freq*t);
    
    % Apply a low-pass filter to the demodulated signal to remove high-frequency components
    demodulated_low = lowpass(demodulated_clean, fc, fs, Steepness=0.99);
    % Store the demodulated, filtered audio for later use
    % audioClean
    audioCleanCell{idx} = demodulated_low; %time domain
    
    audioClean_shift = fftshift(fft(audioReceived))/fs;

    audioCleanCell{idx+length(carrier_freq)} = audioClean_shift;

 % Subplot for the time-domain representation
    figure;
    subplot(2, 1, 1); % Two rows, one column, first subplot
    plot(t, demodulated_low);
    title(sprintf('Time Domain - Carrier %d Hz', carrier_freq));
    xlabel('Time (s)');
    ylabel('Amplitude');

    subplot(2,1,2);
    plot(f, abs(audioClean_shift));
    title('Frequency Domain Representation of Received Audio', carrier_freq);
    xlabel('Frequency (HZ)');
    ylabel('Magnitude');
end

%%
% Optional: Listen to each audioReceived after plotting
% for idx = 1:length(audioReceivedCell)/2
%     sound(audioReceivedCell{idx}, fs);
%     pause(length(audioReceivedCell{idx})/fs + 1);  % Play each sound with a pause
% end
%% 1.5 denoise the signal tone frequency
noise_freq = 2040;



%% helper functions

% function definitions in matlab either need to be in their own file,
% or can be in at the bottom of a script.
% if you want to these functions outside this lab, feel free to 
% move them into their own file, just make sure the filename is the same 
% as the function name, ie timevec.m, freqvec.m


function t= timevec(t0, t0_plus_T, n)
% Creates time vector, where upper limit is non-inclusive
%          t0 <= t < t0_plus_T
%   It is the responsibility of the user to ensure that for the use-case
%   that they want the lower limit included, and the upper-limit 
%   not included.
%
%   Args:
%   t0 = start time
%   t0_plus_T = end time (t0 + T)
%   n = number of samples

    t = linspace(t0, t0_plus_T, n + 1);
    t = t(1:end - 1);
end

%% 
function f=freqvec(fs, n)
% Creates frequency vector suitable for plotting magnitude/phase spectrum
%
%  This function emulates the np.fft.freqvec function from python, but will
%  also make sure that the frequency vector has been shifted correctly, so
%  that the first index is for the lowest frequency, highest index is for
%  highest and that the middle frequency is DC.

%  https://numpy.org/doc/stable/reference/generated/numpy.fft.fftfreq.html
%  
%  The frequency vector will be slightly different if the sequence is of 
%  even length or odd length. The specifics of this has to do with 
%  properties of the Discrete Fourier Transform (DTF), which is what the
%  FFT algorigthm actually computes. Will likely cover
%  this in your Digital signal processing class. The main thing to know is 
%  that this function will create a frequency vector that will ensure the DC
%  Component is at the correct location. For now, can just take our
%  word for it, and know what the function does and what it will return.
%
% 
%  For even length signals, our frequency vector will be of the form,
%        -fs/2 <= f < fs/2
%  For odd length signals, will be,
%        -fs/ 2 < f < fs/2
%
%  Args:
%  fs = sample frequency in Hz
%  n = length of the time vector/number of samples


    % if is an even sequence length, generating the frequency vector
    % is just like doing it for our time vector in the timevec function
    if mod(n, 2) == 0
        f_str = sprintf('Generating freq vec\n [%.2f, %.2f)\n', -fs/2, fs/2);
        disp(f_str);
        % compute the frequency vector
        f = linspace(-fs / 2, fs / 2, n + 1);
        f = f(1:end - 1); 
    % otherwise is of odd length
    else        
        f_str = sprintf('Generating freq vec\n (%.2f, %.2f)\n', ...
            -(n -1)/2 * fs / n, (n -1)/2 * fs / n);
        disp(f_str);
        % compute the frequency vector
        f = linspace(-(n -1)/ 2, (n -1) / 2, n) * fs / n;
    end
end