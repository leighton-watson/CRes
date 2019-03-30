%% EXAMPLE 2 %%
%
% Example script file for CRes (Crater Resonance). Compute the infrasound
% signal in the time and frequency domain for a given crater geometry and
% specified source-time function as a function of crater depth/position of
% lava lake.

clear all;
clc;
cmap = get(gca,'ColorOrder');

% set default plotting options
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',16);

% add file paths
addpath ../source/resonance
addpath ../source/SBPoperators/

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
figure(2); subplot(1,3,1);
plot(shape(:,2), shape(:,1),'k');
set(gca,'YDir','Reverse');
hold on;
plot(-shape(:,2), shape(:,1),'k');
ylim([0 150]);
xlim([-100 100])
xlabel('Radius (m)');
ylabel('Depth (m)');
axis normal

return


%% RESONANCE 1D %%
%
% Iterate over a range of depths

dmin = 20; % min depth
dmax = 140; % max depth
subplot(1,3,1);
ylim([dmin dmax]);
dz = 4; % depth increment
%depths = dmin:dz:dmax; % range of depths
depths = [50 80 120]

depths_plot = [50 80 120]; % plot infrasound signal for select depths
num_depths_plot = length(depths_plot); % number of depth intervals to plot
k = 1; % counter

for i = 1:length(depths)
    
    depth = depths(i);
    str = strcat('Depth=',num2str(depth),'m');
    disp(str); % display string to keep track of iterations while they run
    
    A = resonance1d(shape, depth, freq, Nf, style, order, M); % solve transfer function
    
    transFunc = A.P;
    dp = transFunc .* S(1:N/2+1);
    
    % Compute resonant frequency and quality factor of spectral peak
    [f0(i),Q(i)] = resPeakProps(A.f,dp);
    
    if ismember(depth,depths_plot) % if member depth is member of vector depth plot then plot infrasound signal
    
        figure(3);
        % plot infrasound signal in frequency domain
        subplot(num_depths_plot,2,2*k);
        plot(A.f, abs(dp)./max(abs(dp)),'Color',cmap(k,:));
        xlim([0 3]);
        
        % invert infrasound signal to time domain
        dp_pos = dp(1:N/2+1);
        dp_full = [dp_pos conj(dp_pos(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
        dp_full(N/2+1) = real(dp_full(N/2+1)); % entry at Nyquist must be real
        dp_time = ifft(dp_full,'symmetric')/dt; % inverse fourier transform

        subplot(num_depths_plot,2,(2*k)-1);
        plot(t-tshift, dp_time./max(dp_time),'Color',cmap(k,:));
        xlim([0 10])
        ylim([-1.5 1.1])
               
        k = k+1; % update counter
        
    end
           
end
   
% plot depths_plot
figure(2);
subplot(1,3,1);
hline(depths_plot);
% plot resonant frequency
subplot(1,3,2);
plot(f0,depths,'k');
set(gca,'YDir','Reverse');
xlabel('Frequency (Hz)');

% plot quality factor
subplot(1,3,3);
plot(Q,depths,'k');
set(gca,'YDir','Reverse');
xlabel('Quality Factor');
