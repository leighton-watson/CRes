%% FIG 2 %%
%
% Compute spectrogram for infrasound time series recorded at station EMFO.
% Calculate peak frequency and quality factor in hour long windows. Compute
% average spectra for before and after the onset of the flank eruption.
%
% Written by Leighton Watson
% March 10, 2020
% leightonwatson@stanford.edu // leightonmwatson@gmail.com

clear all; clc;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',1);
set(0,'DefaultAxesFontSize',18);

path(pathdef)
%addpath ../Data/corrected_data
addpath data/
addpath ../source/resonance/

figHand1 = figure(1); clf
set(figHand1,'Position',[100 100 1100 300]);
figHand2 = figure(2); clf
set(figHand2,'Position',[100 100 1100 400]);
figHand3 = figure(3); clf
set(figHand3,'Position',[100 100 1100 300]);


%% load data in to MATLAB %%

disp('Loading infrasound data: Station EMFO');

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
Fs0 = 100; % sampling frequency is 100 Hz
tskip = 10; % number of time steps to skip
Fs = Fs0/tskip; % new sampling frequency
dt = (1/Fs0)*tskip; % time interval
data = data(1:tskip:end); % downsampled data

% generate time vector
k = length(data); % length of data
t = 0:dt:(k-1)*dt; % time vector
t_hrs = t/3600; % time vector in hours

wvs = detrend(data,'constant')*ac_calib;
wvs_filt = wvs;

%% calculate and plot a spectrogram %%

figure(1);
wl = Fs*300; % window length (in seconds)
ol = round(wl*.9); % overlap
nfft = 1024*2;
[S,F,T,P] = spectrogram(wvs,wl,ol,nfft,Fs);
spectrogramPlot = 10*log10(imgaussfilt(P,3));

% normalize spectra at each time
specNorm = zeros(size(P));
for i = 1:length(T)
    i
    specNorm(:,i) = P(:,i)./max(P(:,i));
end

specNormPlot = (imgaussfilt(specNorm,10));
h = surf(T/3600,F,specNormPlot);
view(2); shading interp
ylim([0 3])
xlim([0 72])
h.Parent.XTick = [0 12 24 36 48 60 72];
title('Station EMFO');
xlabel('Time (hr)');
ylabel('Frequency (Hz)');

%% compute entire time spectra

erupt_start = 35; % start of eruption
erupt_end = 38; % end of eruption

% compute spectra before
start_idx = find(t_hrs < erupt_start, 1, 'last');

t_tmp = t_hrs(1:start_idx);
data_tmp = wvs(1:start_idx);

L = length(data_tmp);
NFFT = 2^nextpow2(L);
DataBefore = fft(data_tmp,NFFT)/L;
DataBeforeNorm = abs(DataBefore)./max(abs(DataBefore));
fB = Fs/2*linspace(0,1,NFFT/2+1);
nfilt = 1000;
DataBeforePlot = medfilt1(DataBeforeNorm(1:NFFT/2+1),nfilt);
DataBeforePlot = DataBeforePlot./max(DataBeforePlot);

figure(2); clf;
subplot(1,2,1);
h = plot(fB, DataBeforePlot,'k');
h.Parent.YTick = [0 0.5 1];
xlim([0 3])
hold on;
xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
title('Before eruption onset');


% compute spectra after
stop_idx = find(t_hrs > erupt_end, 1, 'first');

t_tmp = t_hrs(stop_idx:end);
data_tmp = wvs(stop_idx:end);

L = length(data_tmp);
NFFT = 2^nextpow2(L);
DataAfter = fft(data_tmp,NFFT)/L;
DataAfterNorm = abs(DataAfter)./max(abs(DataAfter));
fA = Fs/2*linspace(0,1,NFFT/2+1);
DataAfterPlot = medfilt1(DataAfterNorm(1:NFFT/2+1),nfilt);
DataAfterPlot = DataAfterPlot./max(DataAfterPlot);

subplot(1,2,2);
h = plot(fA, DataAfterPlot,'k');
h.Parent.YTick = [0 0.5 1];
xlim([0 3])
hold on;
xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
title('After eruption onset');

%% short time windows

%%% before %%%
before0 = 15;
before1 = 15.25;

% define bounds of window
idx1 = find(t_hrs > before0, 1, 'first');
idx2 = find(t_hrs < before1, 1, 'last');

% select time and data within the window
t_before = t_hrs(idx1:idx2);
data_before = wvs_filt(idx1:idx2);

% Fourier transform
L = length(data_before);
NFFT2 = 2^nextpow2(L);
Data_freq_tmp = fft(data_before,NFFT2)/L;
f2 = Fs/2*linspace(0,1,NFFT2/2+1);

% save spectra
DataBefore = Data_freq_tmp(1:NFFT2/2+1); % spectra
DataBefore = medfilt1(abs(DataBefore),100);
DataBeforeNorm = DataBefore./max(DataBefore);
subplot(1,2,1);
plot(f2, DataBeforeNorm,'Color',cmap(1,:),'LineWidth',2)
[f0_before, Q_before] = resPeakProps_find(f2, DataBeforeNorm);
legend('Entire period','15 min of high SNR');

%%% after %%%%
after0 = 48 + 6;
after1 = 48 + 6.25;

% define bounds of window
idx1 = find(t_hrs > after0, 1, 'first');
idx2 = find(t_hrs < after1, 1, 'last');

% select time and data within the window
t_after = t_hrs(idx1:idx2);
data_after = wvs_filt(idx1:idx2);

% Fourier transform
L = length(data_after);
NFFT2 = 2^nextpow2(L);
Data_freq_tmp = fft(data_after,NFFT2)/L;
f2 = Fs/2*linspace(0,1,NFFT2/2+1);

% save spectra
DataAfter = Data_freq_tmp(1:NFFT2/2+1); % spectra
DataAfter = medfilt1(abs(DataAfter),100);
DataAfterNorm = DataAfter./max(DataAfter);
subplot(1,2,2);
plot(f2, DataAfterNorm,'Color',cmap(2,:),'LineWidth',2)
[f0_after, Q_after] = resPeakProps_find(f2, DataAfterNorm);
legend('Entire period','15 min of high SNR');

%%% save time series and spectra %%%
EMFO_before.t = t_before;
EMFO_before.p = data_before;
EMFO_before.f = f2;
EMFO_before.spec = DataBeforeNorm;
EMFO_before.f0 = f0_before;
EMFO_before.Q = Q_before;

EMFO_after.t = t_after;
EMFO_after.p = data_after;
EMFO_after.f = f2;
EMFO_after.spec = DataAfterNorm;
EMFO_after.f0 = f0_after;
EMFO_after.Q = Q_after;

save('EMFO_spectra','EMFO_before','EMFO_after')

%% Station EPCN 
disp('Loading infrasound data: Station ECPN');

ac_calib = 0.001/50; % convert from counts to Pa 

disp('December 23');
data1 = load('181223_0000_ECPNo_corrected.mat');
data1 = data1.signal_tot;

% data = data1;
disp('December 24');
data2 = load('181224_0000_ECPNo_corrected.mat');
data2 = data2.signal_tot;

disp('December 25');
data3 = load('181225_0000_ECPNo_corrected.mat');
data3 = data3.signal_tot;

data = [data1; data2; data3];
disp('Finished loading data');

% downsample data to 10 Hz sampling rate
Fs0 = 100; % sampling frequency is 100 Hz
tskip = 10; % number of time steps to skip
Fs = Fs0/tskip; % new sampling frequency
dt = (1/Fs0)*tskip; % time interval
data = data(1:tskip:end); % downsampled data

% generate time vector
k = length(data); % length of data
t = 0:dt:(k-1)*dt; % time vector
t_hrs = t/3600; % time vector in hours

wvs = detrend(data,'constant')*ac_calib;
wvs_filt = wvs;

%% calculate and plot a spectrogram %%

figure(3);
wl = Fs*300; % window length (in seconds)
ol = round(wl*.9); % overlap
nfft = 1024*2;
[S,F,T,P] = spectrogram(wvs,wl,ol,nfft,Fs);
spectrogramPlot = 10*log10(imgaussfilt(P,3));

% normalize spectra at each time
specNorm = zeros(size(P));
for i = 1:length(T)
    i
    specNorm(:,i) = P(:,i)./max(P(:,i));
end

specNormPlot = (imgaussfilt(specNorm,10));
h = surf(T/3600,F,specNormPlot);
view(2); shading interp
ylim([0 3])
xlim([0 72])
h.Parent.XTick = [0 12 24 36 48 60 72];
title('Station ECPN');
xlabel('Time (hr)');
ylabel('Frequency (Hz)');
