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
%addpath ../../../Data/corrected_data
%addpath /Users/lwat054/Documents/Research/Etna_Resonance/Data/corrected_data/

figHand1 = figure(1); clf
set(figHand1,'Position',[100 100 900 400]);

%% load data in to MATLAB %%

disp('Loading infrasound data');

ac_calib = 0.001/50; % convert from counts to Pa 

disp('December 23');
data1 = load('181223_0000_EMFOo_corrected.mat');
data1 = data1.signal_tot;

% data = data1;
disp('December 24');
data2 = load('181224_0000_EMFOo_corrected.mat');
data2 = data2.signal_tot;

disp('December 25');
data3 = load('181225_0000_EMFOo_corrected.mat');
data3 = data3.signal_tot;

data = [data1; data2; data3];
disp('Finished loading data');

% downsample data to 10 Hz sampling rate
Fs0 = 100; % sampling frequency is 100 Hz (samples per second)
tskip = 1; % number of time steps to skip
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
ylim([-11 11])

% plot eruption time window
erupt_start = 35; % start of eruption
erupt_end = 38; % end of eruption
h = rectangle('Position',[erupt_start -11 3 22]);
h.EdgeColor = [1 0 0 0.5];
h.FaceColor = [1 0 0 0.5];

% plot high signal to noise periods that we invert
before0 = 15;
before1 = 15.25;

after0 = 48 + 6;
after1 = 48 + 6.25;

h = rectangle('Position',[before0 -11 0.25 22]);
h.EdgeColor = cmap(1,:);
h.EdgeColor(4) = 0.5;
h.FaceColor = cmap(1,:);
h.FaceColor(4) = 0.5;

h = rectangle('Position',[after0 -11 0.25 22]);
h.EdgeColor = cmap(2,:);
h.EdgeColor(4) = 0.5;
h.FaceColor = cmap(2,:);
h.FaceColor(4) = 0.5;

%% SOUND PRESSURE LEVEL

n = length(wvs_filt); % length of time series
sps = Fs; % samples per second
win_length = sps*5*60; % five minute window
win_overlap = 0.90; % overlap between windows
win_diff = 1 - win_overlap; % difference between windows
win_unique = ceil(win_length*win_diff);
num_win = floor(n/win_unique)-10;
p0 = 20e-6; % reference pressure value

for i = 1:num_win
   idx0(i) = 1 + win_unique*(i-1);
   idx1(i) = idx0(i) + win_length;
   idxc(i) = idx0(i) + ceil(win_length/2);
   tc(i) = t_hrs(idxc(i));
   
   p_tmp = wvs_filt(idx0(i):idx1(i));
   prms(i) = sqrt(sum(p_tmp.^2)/length(p_tmp));
    
   SPL(i) = 20*log10(prms(i)/p0);
    
end

subplot(2,1,2);
h = plot(tc, SPL,'k');
xlim([0 72])
ylim([50 100])
h.Parent.XTick = [0 12 24 36 48 60 72];
h.Parent.XTickLabel = {};
ylabel('Sound Pressure Level (dB re 20 \mu Pa)');



% %% calculate RSAM
% 
% rsam_interval = Fs*5*60; % five minute RSAM interval
% cumulative_function = cumsum(abs(wvs_filt));
% rsam_intervals = rsam_interval:rsam_interval:length(cumulative_function);
% cumulative_function_decimated = cumulative_function(rsam_intervals,:);
% rsam_function = diff(cumulative_function_decimated);
% subplot(2,1,2);
% %plot(rsam_intervals(1:end-1),rsam_function,'k')
% h = plot(t(rsam_intervals(1:end-1))/3600,rsam_function,'k');
% h.Parent.XTick = [0 12 24 36 48 60 72];
% xlabel('time (s)')
% xlim([0 72])
% ylabel('arbitrary units')
% ylim([0 25000])
% 
% plot eruption time window
h = rectangle('Position',[erupt_start 0 3 200]);
h.EdgeColor = [1 0 0 0.5];
h.FaceColor = [1 0 0 0.5];

h = rectangle('Position',[before0 0 0.25 200]);
h.EdgeColor = cmap(1,:);
h.EdgeColor(4) = 0.5;
h.FaceColor = cmap(1,:);
h.FaceColor(4) = 0.5;

h = rectangle('Position',[after0 0 0.25 2000]);
h.EdgeColor = cmap(2,:);
h.EdgeColor(4) = 0.5;
h.FaceColor = cmap(2,:);
h.FaceColor(4) = 0.5;

