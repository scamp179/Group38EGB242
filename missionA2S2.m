%% EGB242 Assignment 2, Section 2 %%
% This file is a template for your MATLAB solution to Section 2.
%
% Before starting to write code, generate your data with the ??? as
% described in the assignment task.

%% Initialise workspace
clear all; close all;

% Begin writing your MATLAB solution below this line.

%% 2.1
% samples = 10^4;
% tvec = linspace(0,20, samples +1);
% tvec(end) = [];
% 
% gm = @(t) 2*(-2 + t + 2*exp(-0.5*t));
% gmvec = gm(tvec);
% 
% figure;
% plot(tvec, gmvec);
% title('DC Motor System Step Response Output');
% xlabel('Time (seconds)');
% ylabel('Voltage Output (\Psi)');

%% 2.1

% Constants
Km = 1;          % Motor constant
alpha = 0.5;     % Motor time constant

% Define Transfer Function G(s)
s = tf('s');
G = Km / (s * (s + alpha));

% Time Vector
samples = 10^4;  % Number of samples
t = linspace(0, 20, samples);  % Time vector from 0 to 20 seconds

% Calculate the Step Response of the System
[y, t_out] = step(G, t);

% Create the Step Input
step_input = ones(size(t_out));  % Step input is 1 for all t

% Convert hexadecimal to normalised RGB triplet

% Plotting
figure;
plot(t_out, y, 'c-', 'LineWidth', 2);          % Step response in cyan
hold on;
plot(t_out, step_input, 'k-', 'LineWidth', 2); % Step input in black line
title('Comparison of Step Input and Step Response');
xlabel('Time (seconds)');
ylabel('System Output (\Psi)');
legend('Step Response', 'Step Input');

%% 2.2

% Define the transfer function for the feedback system
numerator = [1];
denominator = [1, 0.5, 1];
F = tf(numerator, denominator);

% Display poles of the system
poles = pole(F);
disp('Poles of the system:');
disp(poles);

% Time vector from 0 to 20 seconds with 10^4 samples
t = linspace(0, 20, 10^4);

% Calculate the step response of the feedback system
[y, t_out] = step(F, t);

% Plotting the step response and the step input
figure;
plot(t_out, y, 'c-', 'LineWidth', 2);          % Step response in cyan
hold on;
plot(t_out, step_input, 'k-', 'LineWidth', 2); % Step input in black line
title('Comparison of Step Input and Step Response with Feedback');
xlabel('Time (seconds)');
ylabel('System Output');
legend('Step Response', 'Step Input');

%% 2.3
% Constants
zeta = 0.25;  % Damping ratio
omega_n = 1;  % Natural frequency

% Time to Peak
T_p = pi / (omega_n * sqrt(1 - zeta^2));

% Settling Time (2% criterion)
T_s = 4 / (zeta * omega_n);

% Percentage Overshoot
percentage_OS = 100 * exp(-pi * zeta / sqrt(1 - zeta^2));

% Display the results
fprintf('Natural Frequency: %.2f rad/s\n', omega_n);
fprintf('Damping Ratio: %.2f\n', zeta);
fprintf('Time to Peak: %.2f s\n', T_p);
fprintf('Settling Time: %.2f s\n', T_s);
fprintf('Percentage Overshoot: %.2f%%\n'), percentage_OS

%% 2.4
% System Constants
Gm = tf([1], [1, alpha, 0]);  % Motor and load transfer function Gm(s) = 1 / (s*(s + alpha))
Hp = 1;  % Potentiometer gain (Hp(s) = 1 as Kpot = 1)

% Array of gain values to test
K_values = [0.1, 0.2, 0.5, 1, 2];

% Loop through each scenario
for i = 1:2
    figure;  % Create a new figure for each scenario
    for K_gain = K_values
        if i == 1
            % Scenario 1: Vary K_fb with constant K_fwd = 1
            K_fb = K_gain;
            K_fwd = 1;
            legendInfo = sprintf('K_{fb} = %g', K_fb);
            titleInfo = 'Step Response with Varying K_{fb} (K_{fwd} = 1)';
        else
            % Scenario 2: Vary K_fwd with constant K_fb = 1
            K_fwd = K_gain;
            K_fb = 1;
            legendInfo = sprintf('K_{fwd} = %g', K_fwd);
            titleInfo = 'Step Response with Varying K_{fwd} (K_{fb} = 1)';
        end
        
        % Closed-loop transfer function
        F = Gm * K_fwd / (1 + Gm * Hp * K_fb);
        
        % Compute and plot the step response
        [y, t_out] = step(F, t);
        plot(t_out, y, 'DisplayName', legendInfo);
        hold on;  % Hold on to plot multiple lines
    end
    
    % Formatting the plots
    title(titleInfo);
    xlabel('Time (seconds)');
    ylabel('System Output (\Psi)');
    legend show;  % Show legends
end

%% 2.5 

% % Damping ratio and time to peak
% zeta_2 = 0.7;
% T_p_2 = 13;  % Time to peak in seconds
% 
% % Calculate the natural frequency omega_n
% omega_n = pi / (T_p_2 * sqrt(1 - zeta_2^2));
% 
% % Define motor and load transfer function Gm(s)
% Gm = tf([1], [1, alpha, 0]);
% 
% % Forward and feedback gains, assume ideal values for now
% K_fwd = 1;  % This might be adjusted after further analysis
% K_fb = 1;  % This might be adjusted after further analysis
% 
% % Define the transfer function with gains
% F = Gm * K_fwd / (1 + Gm * 1 * K_fb);  % Using Kpot = 1 for H_p(s)
% 
% % Calculate and set the desired system parameters
% desired_omega_n = omega_n;
% desired_zeta = zeta_2;
% 
% % Store the transfer function object for the camera control system
% cameraTF = F;
% 
% % Display the transfer function
% disp(cameraTF);

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

%% helper functions

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
