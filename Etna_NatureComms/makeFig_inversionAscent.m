%% MAKE FIG INVERSE DEPTH %%
%
% Make figure showing the spectrogram for one infrasound station, peak
% frequency of data as a function of time, and inverted depth as a function
% of time.

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
cmap = get(gca,'ColorOrder');

load inversionAscent.mat

figHand = figure(1); clf; 
set(figHand,'Position',[100 100 900 500]);

% ascent velocity
subplot(2,1,1);
h = plot(time,vel,'.');
xlim([-2 24]);
ylim([-20 40])
grid on;
ylabel('Velocity (m/s)');
xlabel('Time on 20 February 2021 (hours)');
set(h(1).Parent,'XTick',[-2:2:24]);
vline(19); vline(20); vline(22);


% volumetric flux
subplot(2,1,2);
h = plot(time,vol,'.');
xlim([-2 24]);
ylim([-5000 10000]);
grid on;
xlabel('Time on 20 February 2021 (hours)');
ylabel('Volumetric');
set(h(1).Parent,'XTick',[-2:2:24]);
vline(19); vline(20); vline(22);