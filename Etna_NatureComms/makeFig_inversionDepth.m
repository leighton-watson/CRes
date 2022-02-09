%% MAKE FIG INVERSE DEPTH %%
%
% Make figure showing the spectrogram for one infrasound station, peak
% frequency of data as a function of time, and inverted depth as a function
% of time.

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
cmap = get(gca,'ColorOrder');

%% spectra for station ECPN
load specECPN.mat;
figHand = figure(1); clf; 
set(figHand,'Position',[100 100 900 300]);
h = surf(T/3600-24,F,Pnorm);
set(h.Parent,'XTick',[-2:2:24]);
view(2); shading interp
colormap parula
ylim([0 5])
ylabel('Frequency (Hz)');
xlabel('Time on 20 February 2021 (hours)');
xlim([-2 24])
hold on;

%% peak frequency of data as a function of time
load etna2021_data_peak_freq.mat
figHand = figure(2); clf;
set(figHand,'Position',[100 100 900 300]);
for i = 1:5
    h = plot(t-24,f0smooth(i,:),'.','MarkerSize',10);
    hold on;
end
xlim([-2 24]);
ylim([0 3.5]); grid on;
xlabel('Time on 20 February 2021 (hours)');
ylabel('Frequency (Hz)');
set(h.Parent,'XTick',[-2:2:24]);
vline(19); % weak strombolian
vline(20); % strong strombolian
vline(22); % lava fountaining

%% inverted depth as a function of time

load inversionDepth_200C.mat
figHand = figure(3); clf;
set(figHand,'Position',[100 100 900 300]);
for i = 1:5
    h = plot(out1.t-24, out1.D(i,:),'.','MarkerSize',10);
    hold on;
end
xlim([-2 24]);
ylim([50 250]); grid on;
set(gca,'YDir','Reverse');
xlabel('Time on 20 February 2021 (hours)');
ylabel('Frequency (Hz)');
set(h.Parent,'XTick',[-2:2:24]);
vline(19); % weak strombolian
vline(20); % strong strombolian
vline(22); % lava fountaining