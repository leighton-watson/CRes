%% EXAMPLE 1 %%
%
% Example script file for CRes (Crater Resonance). Compute the infrasound
% signal in the time and frequency domain for a given crater geometry and
% specified source-time function.

clear all;
clc;
cmap = get(gca,'ColorOrder');

% set default plotting options
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',16);

% add file paths
addpath ../source/resonance
addpath ../source/SBPoperators

%% RESONANCE 1D PARAMETERS %%
%
% Specify parameters for solver

T = 50; % total time [s]
N = 500; % number of grid points in time (formulas below assume even N)
dt = T/N; % time step
t = [0:N-1]*dt; % vector of time grid points [s]
tshift = T/4;

Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples

M = problemParameters(); % property values are defined in resonance/problemParameters.m
style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method

%% SOURCE FUNCTION %%
%
% Compute source function (volumetric flow rate [m^3/s]) at base of crater
resParams = [T, N];
srcStyle = 'Gauss'; % choose from "Gauss" or "Brune"
A = 1; % amplitude of source
L = 0.2; % width of source 
[S,f,s,ts] = sourceFunction(A,L,srcStyle, resParams);

figure(1); clf;
subplot(1,2,1);
plot(ts-tshift+1,s);
xlabel('Time (s)');
ylabel('s(t)');
xlim([0 5])

subplot(1,2,2);
plot(f(1:N/2+1),abs(S(1:N/2+1))./max(abs(S(1:N/2+1))));
xlabel('Frequency (Hz)');
ylabel('s(\omega)');
xlim([0 3])

%% GEOMETRY %%
%
% Load crater geometry from .mat file

load('Johnson2018'); % load geometry
shape = flipud(geometry); % reformat geometry matrix to be read by solver
figure(2);
plot(shape(:,2), shape(:,1),'k');
set(gca,'YDir','Reverse');
hold on;
plot(-shape(:,2), shape(:,1),'k');
axis equal
ylim([0 150]);
xlabel('Radius (m)');
ylabel('Depth (m)');

%% RESONANCE 1D %%
depth = 120; % depth of crater bottom
A = resonance1d(shape, depth, freq, Nf, style, order, M);

transFunc = A.P;
dp = transFunc .* S(1:N/2+1);

figure(3); clf;
subplot(1,2,1);
plot(A.f(1:N/2+1), abs(dp(1:N/2+1))./max(abs(dp(1:N/2+1))));
xlim([0 3]);
xlabel('Frequency (Hz)');
ylabel('\Delta p(\omega,r)');

% invert infrasound signal to the time domain
dp_pos = dp(1:N/2+1);
dp_full = [dp_pos conj(dp_pos(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
dp_full(N/2+1) = real(dp_full(N/2+1)); % entry at Nyquist must be real
dp_time = ifft(dp_full,'symmetric')/dt; % inverse fourier transform

subplot(1,2,2);
plot(t-tshift, dp_time);
xlim([0 15])
xlabel('Time (s)');
ylabel('\Delta p(t,r)');

%% PROPERTIES OF SPECTRAL PEAK
%
% Compute resonant frequency and quality factor of spectral peak
[f0,Q] = resPeakProps(A.f,dp);
