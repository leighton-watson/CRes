%% MAKE FIG %%
%
% Calculate synthetic infrasound signal for a pipe for different gas
% compositions (R values)

clear all;
clc;
close all;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',16);
cmap = get(gca,'ColorOrder');

% add file paths
addpath ../source/resonance
addpath ../source/SBPoperators
addpath Functions/

figHand = figure(1); clf;
set(figHand,'Position',[100 100 1000 280]);

% crater properties
radius = 100; % crater radius
depth = 300; % crater depth
z = [depth:-1:0]'; % depth vector
shape = [z radius*ones(size(z))]; % crater geometry matrix

% parameters for numerical model
T = 100; % total time [s]
N = 2000; % number of grid points in time (formulas below assume even N)
dt = T/N; % time step
t = [0:N-1]*dt; % vector of time grid points [s]
Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples

% problem parameters % save problem parameters in to the structure M
TA = 0; % atmospheric temperature
TC = 100; % crater temperature
RA = 287; % specific gas constant for atmosphere (dry air)
rhoA = 1; % atmospheric density

% source
amp = 1; % amplitude
L = 0.3; % width
srcStyle = 'Gauss'; % style
resParams = [T N]; % resonance parameters
[S, f, s, ts] = sourceFunction(amp, L, srcStyle, resParams); % compute source

% gas compositions
R_air = 287;
R_CO2 = 189;
R_SO2 = 130;

gamma_air = 1.4;
gamma_CO2 = 1.289;
gamma_SO2 = 1.29;

percent_air = 60;
percent_CO2 = 20;
percent_SO2 = 20;

R_mix = percent_air/100*R_air + percent_CO2/100*R_CO2 + percent_SO2/100*R_SO2;
gamma_mix = percent_air/100*gamma_air + percent_CO2/100*gamma_CO2 + percent_SO2/100*gamma_SO2;

R = [R_air R_mix];
gamma = [gamma_air gamma_mix];

LineStyle = {'-','-'};
alpha = [1 0.9];
LineWidth = [2 3];


for i = 1:length(R)
    
    RC = R(i);
    M = problemParameters_GasComp(gamma(i),TA,TC,RA,RC,rhoA);
        
    % solver
    A = resonance1d(shape, depth, freq, Nf, [], [], M);
    
    B = A.pOutlet(1:N/2+1); % field to plot and compute properties for
    B = B.*S(1:N/2+1);
    C.f = A.f;
    C.P = B;
    
    % resonant frequency and quality factor
    [f0(i) Q(i)] = resPeakProps(C, 'sim');
    
    BF = [B conj(B(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
    BF(N/2+1) = real(BF(N/2+1)); % entry at Nyquist must be real
    B_time = ifft(BF,'symmetric')/dt; % inverse Fourier transform to obtain infrasound signal in time domain
    
    % plotting spectra
    subplot(1,3,3);
    if i == 1
        spectra_max = max(abs(B));
    end
    %h = plot(A.f, abs(B)./max(abs(B)),'LineStyle',LineStyle{i},...
    %    'LineWidth',LineWidth(i));
    h1 = plot(A.f, abs(B)./spectra_max);
    h1.Color(4) = alpha(i);
    xlabel('Frequency (Hz)');
    hold on;
    set(h1.Parent,'YTick',[]);
    xlim([0 1.5])
    ylim([0 1.2])
    ylabel('$\Delta p(\omega,r)$','interpreter','latex')
    
    subplot(1,3,[1 2]);
    %h = plot(t-T/4+5, B_time./max(abs(B_time)),'LineStyle',LineStyle{i},...
    %    'LineWidth',LineWidth(i));
    if i == 1
        p_max = max(abs(B_time));
    end
    h2 = plot(t-23+2.75, B_time./p_max);
    h2.Color(4) = alpha(i);
    xlim([0 50])
    ylim([-0.9 1.2])
    hold on;
    set(h2.Parent,'YTick',[]);
    xlabel('Time (s)');
    ylabel('$\Delta p(t,r)$','interpreter','latex')
end
