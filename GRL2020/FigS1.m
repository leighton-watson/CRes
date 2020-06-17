%% FIG S1 %%
%
% Display unfiltered and filtered infrasound spectra recorded at station
% EMFO for before and after the onset of the flank eruption.
%
% Written by Leighton Watson
% March 11, 2020
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
Fs0 = 100; % sampling frequency is 100 Hz
tskip = 10; % number of time steps to skip
Fs = Fs0/tskip; % new sampling frequency
dt = (1/Fs0)*tskip; % time interval
data = data(1:tskip:end); % downsampled data

% generate time vector
k = length(data); % length of data
t = 0:dt:(k-1)*dt; % time vector
t_hrs = t/3600; % time vector in hours

% eruption time window
erupt_start = 35; % start of eruption
erupt_end = 38; % end of eruption

%% filter data %%

wvs = detrend(data,'constant')*ac_calib;
[B,A] = butter(4,0.25/(Fs/2),'high');
wvs_filt = filter(B,A,wvs);

%% calculate spectra in windows (filtered and unfiltered) %%

number_of_hours = 72; % total number of hours of data
number_of_mins = number_of_hours*60; % total number of minutes of data
window_length = 5; % length of window (in minutes) to compute Fourier transform over
number_windows = number_of_mins/window_length; % number of windows to compute Fourier transform over

f0 = zeros(number_windows,1); % peak frequency
Q = zeros(number_windows,1); % quality factor

start_idx = find(t_hrs < erupt_start, 1, 'last');
stop_idx = find(t_hrs > erupt_end, 1, 'first');

disp('Computing Fourier transforms and generating spectrogram');

for i = 1:number_windows
    
    disp(strcat('Percentage complete:',num2str(i/number_windows*100),'%'));
    
    % define bounds of window
    idx1 = find(t_hrs > (i-1)*number_of_hours/number_windows, 1, 'first');
    idx2 = find(t_hrs < i*number_of_hours/number_windows, 1, 'last');
    
    % select time window
    t_tmp = t_hrs(idx1:idx2);
    
    
    % select data within the window
    data_tmp = wvs(idx1:idx2); % unfiltered
    data_tmp_filt = wvs_filt(idx1:idx2); % filtered
    
    % Fourier transform
    L = length(data_tmp);
    NFFT = 2^nextpow2(L);
    DataF_tmp = fft(data_tmp,NFFT)/L; % unfiltered
    DataF_tmp_filt = fft(data_tmp_filt,NFFT)/L; % filtered
    f = Fs/2*linspace(0,1,NFFT/2+1);
    
    % compute peak frequency and quality factor
    [f0(i), Q(i)] = resPeakProps(f, abs(DataF_tmp(1:NFFT/2+1)));
    [f0_filt(i), Q_filt(i)] = resPeakProps(f, abs(DataF_tmp_filt(1:NFFT/2+1)));
    
    % save spectra as before or after eruption
    if max(t_tmp) < t_hrs(start_idx) % before eruption
        
        % unfiltered
        DataB(:,i) = DataF_tmp(1:NFFT/2+1); % spectra
        DataB_norm(:,i) = abs(DataF_tmp(1:NFFT/2+1))....
            /max(DataF_tmp(1:NFFT/2+1)); % normalized spectra
        
        % filtered
        DataB_filt(:,i) = DataF_tmp_filt(1:NFFT/2+1); % spectra
        DataB_filt_norm(:,i) = abs(DataF_tmp_filt(1:NFFT/2+1))....
            /max(DataF_tmp_filt(1:NFFT/2+1)); % normalized spectra
        
    elseif min(t_tmp) > (t_hrs(end) - t_hrs(stop_idx)) % after eruption
        
        % unfiltered
        DataA(:,i) = DataF_tmp(1:NFFT/2+1); % spectra
        DataA_norm(:,i) = abs(DataF_tmp(1:NFFT/2+1))....
            /max(DataF_tmp(1:NFFT/2+1)); % normalized spectra
        
        % filtered
        DataA_filt(:,i) = DataF_tmp_filt(1:NFFT/2+1); % spectra
        DataA_filt_norm(:,i) = abs(DataF_tmp_filt(1:NFFT/2+1))....
            /max(DataF_tmp_filt(1:NFFT/2+1)); % normalized spectra
        
    end
       
end


%% display average spectra before and after flank eruption onset %%

%%% before eruption %%%
DataB_stack = mean(abs(DataB_norm),2);
DataB_stack = DataB_stack./max(DataB_stack);

figure(1); clf;
subplot(1,2,1);
plot(f, DataB_stack','Color',cmap(2,:));
xlabel('Frequency (Hz)');
xlim([0 3]);
ylabel('Normalized Amplitude');
hold on;
title('Before Eruption');

DataB_filt_stack = mean(abs(DataB_filt_norm),2);
DataB_filt_stack = DataB_filt_stack./max(DataB_filt_stack);
plot(f, DataB_filt_stack,'k');

%%% after eruption %%%
DataA_stack = mean(abs(DataA_norm),2);
DataA_stack = DataA_stack./max(DataA_stack);

subplot(1,2,2);
plot(f, DataA_stack,'Color',cmap(2,:));
xlabel('Frequency (Hz)');
xlim([0 3]);
ylabel('Normalized Amplitude');
hold on;
title('After Eruption');

DataA_filt_stack = mean(abs(DataA_filt_norm),2);
DataA_filt_stack = DataA_filt_stack./max(DataA_filt_stack);
plot(f, DataA_filt_stack,'k');

%%% compute peak frequency and quality factor from average spectra
[f0B_avgSpec, QB_avgSpec] = resPeakProps(f, DataB_stack);
[f0B_filt_avgSpec, QB_filt_avgSpec] = resPeakProps(f, DataB_filt_stack);

[f0A_avgSpec, QA_avgSpec] = resPeakProps(f, DataA_stack);
[f0A_filt_avgSpec, QA_filt_avgSpec] = resPeakProps(f, DataA_filt_stack);




