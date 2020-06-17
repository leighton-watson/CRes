%% FIG 1 %%
%
% Display infrasound time series recorded at station EMFO. Compute real
% time seismic amplitude measurement.
%
% Written by Leighton Watson
% March 10, 2020
% leightonwatson@stanford.edu // leightonmwatson@gmail.com

clear all; clc;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',1);
set(0,'DefaultAxesFontSize',18);

path(pathdef)
addpath data/
addpath ../source/resonance/
addpath ../source/SBPoperators/

figHand1 = figure(1); clf
set(figHand1,'Position',[100 100 900 400]);

%% load data in to MATLAB %%

disp('Loading infrasound data');

ac_calib = 0.001/50; % convert from counts to Pa 

disp('December 23');
data1 = load('181223_EMFOo.mat');
data1 = data1.signalz_tot;

% data = data1;
disp('December 24');
data2 = load('181224_EMFOo.mat');
data2 = data2.signalz_tot;

disp('December 25');
data3 = load('181225_EMFOo.mat');
data3 = data3.signalz_tot;

data = [data1; data2; data3];
disp('Finished loading data');

% downsample data to 10 Hz sampling rate
Fs0 = 100; % sampling frequency is 100 Hz (samples per second)
tskip = 10; % number of time steps to skip
Fs = Fs0/tskip; % new sampling frequency (samples per second)
dt = (1/Fs0)*tskip; % time interval
data = data(1:tskip:end); % downsampled data

% generate time vector
k = length(data); % length of data
t = 0:dt:(k-1)*dt; % time vector
t_hrs = t/3600; % time vector in hours

%% filter and plot data time series 

figure(1); clf;
subplot(2,1,1);
wvs = detrend(data,'constant')*ac_calib;
[B,A] = butter(4,0.25/(Fs/2),'high');
wvs_filt = filter(B,A,wvs);
h = plot(t_hrs, wvs_filt,'k');
ylabel('Amplitude (Pa)');
xlim([min(t_hrs) max(t_hrs)]);
h.Parent.XTick = [0 12 24 36 48 60 72];
h.Parent.XTickLabel = {};
ylim([-8 8])

% plot eruption time window
erupt_start = 35; % start of eruption
erupt_end = 38; % end of eruption
h = rectangle('Position',[erupt_start -8 3 16]);
h.EdgeColor = [1 0 0 0.5];
h.FaceColor = [1 0 0 0.5];

%% calculate RSAM

rsam_interval = Fs*5*60; % five minute RSAM interval
cumulative_function = cumsum(abs(wvs_filt));
rsam_intervals = rsam_interval:rsam_interval:length(cumulative_function);
cumulative_function_decimated = cumulative_function(rsam_intervals,:);
rsam_function = diff(cumulative_function_decimated);
subplot(2,1,2);
%plot(rsam_intervals(1:end-1),rsam_function,'k')
h = plot(t(rsam_intervals(1:end-1))/3600,rsam_function,'k');
h.Parent.XTick = [0 12 24 36 48 60 72];
xlabel('time (s)')
xlim([0 72])
ylabel('arbitrary units')
ylim([0 2200]);

% plot eruption time window
h = rectangle('Position',[erupt_start 0 3 2200]);
h.EdgeColor = [1 0 0 0.5];
h.FaceColor = [1 0 0 0.5];

