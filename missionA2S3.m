%%%% EGB242 Assignment 2, Section 3 %%
% This file is a template for your MATLAB solution to Section 3.
%
% Before starting to write code, generate your data with the ??? as
% described in the assignment task.

%% Initialise workspace
clear all; close all;
load DataA2 imagesReceived;

% Begin writing your MATLAB solution below this line.
%% 3.1
im1D = imagesReceived(1,:);
numRows = 480;
numCols = 640;

% Converting a received pixel stream to an image matrix
im2D = reshape(im1D, numRows, numCols);

% Displaying an image in a figure
figure;
imshow(im2D ) ;

% Saving an image matrix as an image file (to include in report)
imwrite(im2D , 'pic1.png' ) ;

%% 3.2

fs = 1000;
T = length(im1D) / fs;
N = length(im1D);
t = timevec(0,T,N);
f = freqvec(fs, N);

shift_im1D = fftshift(fft(im1D)) / fs;

figure;

subplot(2, 1, 1);
plot(t, im1D);
title('Time Domain Representation of Image Signal');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 1, 2);
plot(f, abs(shift_im1D));
title('Frequency Domain Representation of Image Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');

%% 3.3
R = 820;
C = 1e-6;
s = tf('s');

num = [1 / (R*C)^2];
den = [1, 2/(R*C), 1 / (R*C)^2];
transferFunction = (1 / (R*C)^2)/(s^2 + 2*s/(R*C) + 1/(R*C)^2);
H = tf(num, den);
ltiview(H)

num2 = [1, 0, 0];
transferFunction2 = s^2/(s^2 + 2*s/(R*C) + 1/(R*C)^2);
H2 = tf(num2, den);
ltiview(H2)

%Use Active Filter 1 to remove noise
%% 3.4
G = @(s) (1 / (R*C)^2)/(s^2 + 2*s/(R*C) + 1/(R*C)^2);
syms t v s
V = v*exp(-s*t);
Vs = int(V)
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