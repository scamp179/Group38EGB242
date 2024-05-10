%% EGB242 Assignment 2, Section 2 %%
% This file is a template for your MATLAB solution to Section 2.
%
% Before starting to write code, generate your data with the ??? as
% described in the assignment task.

%% Initialise workspace
clear all; close all;

% Begin writing your MATLAB solution below this line.
samples = 10^4;
tvec = linspace(0,20, samples +1);
tvec(end) = [];

gm = @(t) 2*(-2 + t + 2*exp(-0.5*t));
gmvec = gm(tvec);

figure;
plot(tvec, gmvec);
title('DC Motor System Step Response Output');
xlabel('Time (seconds)');
ylabel('Voltage Output (\Psi)');

