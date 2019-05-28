%% makeFig2

clear all;
close all;
clc;

% set default options
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
figHand1 = figure(1); clf;
set(figHand1, 'Position', [100 100 1600 700]);
cmap = get(gca,'ColorOrder');
alpha = 0.7; % plot transparency
alpha2 = 0.9;

% add file paths
addpath ../source/resonance
addpath ../source/SBPoperators

%% Resonance 1D Parameters %%
T = 100;
N = 2000;
%T = 200; % total time [s]
%N = 4000; % number of grid points in time (formulas below assume even N)
dt = T/N; % time step
t = [0:N-1]*dt; % vector of time grid points [s]

Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples

M = problemParameters(); % property values are defined in resonance/problemParameters.m
style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method

%% Crater Geometry %%

depth_vector = [1000:-1:0]'; 
radius_vector = 100*ones(size(depth_vector));
shape = [depth_vector radius_vector];
depth = 200;

%% Source Function %%

numSrc = 4; % number of source functions to consider
amp = [1 1 1 -1]; % amplitude
L = [0.5 0.2 0.2 0.3]; % width
srcStyle = {'Gauss','Gauss','Brune','Gauss'}; % style
resParams = [T N]; % resonance parameters
tshift = [0.25 0 -0.25 0.25];
%% Resonance 1D %%

for i = 1:numSrc
    
    % compute source
    [S, f, s, ts] = sourceFunction(amp(i), L(i), srcStyle{i}, resParams);
    
    % plot source in time domain
    subplot(numSrc,5,5*(i)-4);
    h = plot(ts-T/4,s./max(abs(s)));
    h.Color(4) = alpha2;
    set(h.Parent,'YTick',[]);
    xlim([-2 3])
    xlabel('Time (s)');
    ylabel('$s(t)$','interpreter','latex')
    
    % plot source in frequency domain
    %subplot(numSrc+1,4,4*(i+1));
    subplot(numSrc,5,5*(i)-3);
    h = plot(f(1:N/2+1), abs(S(1:N/2+1))./max(abs(S(1:N/2+1))));
    h.Color(4) = alpha2;
    set(h.Parent,'YTick',[]);
    xlim([0 2]);
    xlabel('Frequency (Hz)');
    ylabel('$s(\omega)$','interpreter','latex')
    
    % compute Green's function
    A = resonance1d(shape, depth, freq, Nf, style, order, M);
    
    % convolve Green's function with source
    pres = A.P(1:N/2+1) .* S(1:N/2+1);
    
    %%% plot Green's function %%%
    subplot(numSrc,5,5*i);
    greens_func = A.P(1:N/2+1);
    h = plot(A.f(1:N/2+1),abs(greens_func)./max(abs(greens_func)),...
        'Color',cmap(2,:),'LineStyle','-');
    h.Color(4) = alpha;
    set(h.Parent,'YTick',[]);
    xlim([0 2]);
    ylim([0 1]);
    hold on;
    xlabel('Frequency (Hz)');
    ylabel('$\Delta p(\omega,r)$','interpreter','latex')

    greens_funcF = [greens_func conj(greens_func(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
    greens_funcF(N/2+1) = real(greens_funcF(N/2+1)); % entry at Nyquist must be real
    greens_func_time = ifft(greens_funcF,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain
    
    subplot(numSrc,5,[5*i-2 5*i-1]);
    h = plot(t, amp(i)*greens_func_time./max(abs(greens_func_time)),...
        'Color',cmap(2,:),'LineStyle','-');
    h.Color(4) = alpha;
    set(h.Parent,'YTick',[]);
    xlim([0 30]);
    ylim([-1 1]);
    hold on;
    xlabel('Time (s)');
    ylabel('$\Delta p(t,r)$','interpreter','latex')

    
    % plot synthetic infrasound signal in frequency domain
    subplot(numSrc,5,5*i);
    h = plot(A.f(1:N/2+1), abs(pres)./max(abs(pres)),'Color',cmap(1,:));
    h.Color(4) = alpha2;
    set(h.Parent,'YTick',[]);
    xlim([0 2]);
    
    % invert synthetic infrasound signal to time domain
    presF = [pres conj(pres(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
    presF(N/2+1) = real(presF(N/2+1)); % entry at Nyquist must be real
    presTime = ifft(presF,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain

    % plot synthetic infrasound signal in time domain
    subplot(numSrc,5,[5*i-2 5*i-1]);
    h = plot(t-T/4+tshift(i), presTime./max(abs(presTime)),'Color',cmap(1,:));
    h.Color(4) = alpha2;
    set(h.Parent,'YTick',[]);
    xlim([0 30])
%     
%     % label plots
%     if i == numSrc
%         %subplot(numSrc+1,4,4*(i+1)-3);
%         subplot(numSrc+1,6,[6*(i+1)-5 6*(i+1)-4]);
%         xlabel('Time (s)');
%         
%         %subplot(numSrc+1,4,4*(i+1)-2);
%         subplot(numSrc+1,6,[6*(i+1)-3 6*(i+1)-2]);
%         xlabel('Frequency (Hz)');
%         
%         %subplot(numSrc+1,4,4*(i+1)-1);
%         subplot(numSrc+1,6,6*(i+1)-1);
%         xlabel('Time (s)');
%         
%         %subplot(numSrc+1,4,4*(i+1));
%         subplot(numSrc+1,6,6*(i+1));
%         xlabel('Frequency (Hz)');
%         
%     end
%     

    
end




