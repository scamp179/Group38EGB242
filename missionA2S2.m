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
plot(t_out, y);          % Step response in cyan
hold on;
plot(t_out, step_input); % Step input in black line
title('Comparison of Step Input and Step Response');
xlabel('Time (seconds)');
ylabel('System Output (\Psi)');
legend({'Step Response', 'Step Input'}, 'Location', 'northeast');


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
plot(t_out, y);         
hold on;
plot(t_out, step_input); 
title('Comparison of Step Input and Step Response with Feedback');
xlabel('Time (seconds)');
ylabel('System Output');
legend({'Step Response', 'Step Input'}, 'Location', 'northeast');

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

Gm = tf([1], [1, alpha, 0]); % Transfer function Gm(s) = 1 / (s*(s + alpha))
Hp = 1; % Potentiometer gain

K_fb_values = [0.1, 0.2, 0.5, 1, 2]; % Feedback gain values
K_fwd_values = [0.1, 0.2, 0.5, 1, 2]; % Forward gain values



% Scenario 1: Varying K_fb with K_fwd = 1
K_fwd = 1; % Constant forward gain
figure;
for K_fb = K_fb_values
    F = Gm * K_fwd / (1 + Gm * Hp * K_fb); % Transfer function with feedback
    [y, t_out] = step(F, t);
    plot(t_out, y, 'DisplayName', sprintf('K_{fb} = %.1f', K_fb));
    hold on;
end
title('Step Response with Varying K_{fb}, K_{fwd} = 1');
xlabel('Time (seconds)');
ylabel('System Output (\Psi)');
legend('show', 'Location', 'best');


% Scenario 2: Varying K_fwd with K_fb = 1
K_fb = 1; % Constant feedback gain
figure;
for K_fwd = K_fwd_values
    F = Gm * K_fwd / (1 + Gm * Hp * K_fb); % Transfer function with forward gain
    [y, t_out] = step(F, t);
    plot(t_out, y, 'DisplayName', sprintf('K_{fwd} = %.1f', K_fwd));
    hold on;
end
title('Step Response with Varying K_{fwd}, K_{fb} = 1');
xlabel('Time (seconds)');
ylabel('System Output (\Psi)');
legend('show');


%% 2.5 


% Specifications
zeta_2 = 0.7;           % Chosen damping ratio for controlled response
T_p_2 = 13;             % Time to peak in seconds
omega_n_2 = pi / (T_p_2 * sqrt(1 - zeta_2^2));  % Calculate natural frequency

% Closed-Loop Transfer Function
% Redefine Gm(s) to include omega_n and damping ratio zeta
Gm_adjusted = tf([omega_n_2^2], [1, 2*zeta_2*omega_n_2, omega_n_2^2]);
F = Gm_adjusted * K_fwd / (1 + Gm_adjusted * Hp * K_fb);  % Adjusted system

% Simulation
[y, t_out] = step(F, t);   % Step response
% 
% Plotting
figure;
plot(t_out, y);
title('Controlled Step Response of Camera System');
xlabel('Time (seconds)');
ylabel('Camera Position (Radians)');

% Storing the Transfer Function
cameraTF = F;  % Save transfer function for future use
disp('Configured Camera Transfer Function:');
disp(cameraTF);
save('cameraTF.mat', 'cameraTF');

%% 2.6

% Constants
deg_to_rad = pi / 180;  % Convert degrees to radians
start_angle_deg = 30;   % Start angle in degrees
end_angle_deg = 210;    % End angle in degrees

% Convert angles to voltage (assuming linear relationship: 0V at 0°, 1V at 360°)
start_voltage = start_angle_deg / 360;
end_voltage = end_angle_deg / 360;

% Ensure cameraTF is configured and loaded
load cameraTF;  % This line assumes cameraTF is saved and needs to be loaded

% Simulate the camera pan
[startIm, finalIm] = cameraPan(start_voltage, end_voltage, cameraTF);

% Display results
figure;
subplot(1,2,1);
imshow(startIm);
title('Starting Image at 30°');

subplot(1,2,2);
imshow(finalIm);
title('Final Image at 210°');

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
