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
addpath data/
addpath ../source/resonance/
addpath ../source/SBPoperators/

figHand1 = figure(1); clf
set(figHand1,'Position',[100 100 900 1200]);

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

%% filter time series %%

wvs = detrend(data,'constant')*ac_calib;
[B,A] = butter(4,0.25/(Fs/2),'high');
wvs_filt = filter(B,A,wvs);

%% calculate and plot a spectrogram %%

subplot(4,2,[1 2]);
wl = Fs*300; % window length (in seconds)
ol = round(wl*.9); % overlap
nfft = 1024*2;
[S,F,T,P] = spectrogram(wvs_filt,wl,ol,nfft,Fs);
spectrogramPlot = 10*log10(imgaussfilt(P,3));
h = imagesc(T/3600,F,spectrogramPlot);
h.Parent.XTick = [0 12 24 36 48 60 72];
h.Parent.XTickLabel = {};
set(gca,'ydir','normal')
%h = colorbar;
%ylabel(h, 'Spectral Amplitude (dB)')
caxis([-40 0])
ylabel('Frequency (Hz)')
ylim([0 3])
colormap(jet)
hold on;

% plot eruption time window
erupt_start = 35; % start of eruption
erupt_end = 38; % end of eruption
h = patch([erupt_start erupt_end erupt_end erupt_start],[0 0 3 3],[1 1 1 1]*1e3,'r');
h.EdgeColor = [1 0 0];
h.EdgeAlpha = 0.4;
h.FaceColor = [1 0 0];
h.FaceAlpha = 0.4;

%% compute spectra in moving windows %%

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
    
    % select time and data within the window
    t_tmp = t_hrs(idx1:idx2);
    data_tmp = wvs_filt(idx1:idx2);
    
    % Fourier transform
    L = length(data_tmp);
    NFFT = 2^nextpow2(L);
    Data_freq_tmp = fft(data_tmp,NFFT)/L;
    f = Fs/2*linspace(0,1,NFFT/2+1);
    
    % save spectra
    Data_save(:,i) = Data_freq_tmp(1:NFFT/2+1); % spectra
    Data_norm_save(:,i) = abs(Data_freq_tmp(1:NFFT/2+1))./max(Data_freq_tmp(1:NFFT/2+1)); % normalized spectra
        
    % compute peak frequency and quality factor
    [f0(i), Q(i)] = resPeakProps(f, abs(Data_freq_tmp(1:NFFT/2+1)));
    
    % check if time is before or after eruption
    if max(t_tmp) < t_hrs(start_idx) % before eruption
        Data_before_save(:,i) = Data_freq_tmp(1:NFFT/2+1); % spectra
        Data_before_norm_save(:,i) = abs(Data_freq_tmp(1:NFFT/2+1))....
            /max(Data_freq_tmp(1:NFFT/2+1)); % normalized spectra
      
    elseif min(t_tmp) > (t_hrs(end) - t_hrs(stop_idx)) % after eruption 
        Data_after_save(:,i) = Data_freq_tmp(1:NFFT/2+1); % spectra
        Data_after_norm_save(:,i) = abs(Data_freq_tmp(1:NFFT/2+1))....
            /max(Data_freq_tmp(1:NFFT/2+1)); % normalized spectra
    end
       
       
end

%% peak frequency and quality factor %%

time = ((1:number_windows)*number_of_hours/number_windows)'; % time vector
filt_length = 60/window_length;
f0_filt = medfilt1(f0,filt_length);
Q_filt = medfilt1(Q,filt_length);

% peak frequency
figure(1); subplot(4,2,[3 4]); hold on;
plot([0 erupt_start],[0.60 0.60],'LineWidth',3,'Color',cmap(2,:),'LineStyle','--');
plot([erupt_end max(time)],[0.36 0.36],'LineWidth',3,'Color',cmap(3,:),'LineStyle','--');
h = plot(time(3:end), f0_filt(3:end),'Color','k'); 
h.LineWidth = 2;
ylabel('Peak Frequency (Hz)');
xlim([min(t_hrs) max(t_hrs)]);
h.Parent.XTick = [0 12 24 36 48 60 72];
h.Parent.XTickLabel = {};
ylim([0 1.2]);
box on;

% plot eruption time window
h = rectangle('Position',[erupt_start 0 3 1.2]);
h.EdgeColor = [1 0 0 0.5];
h.FaceColor = [1 0 0 0.5];

% quality factor
subplot(4,2,[5 6]); hold on;
plot([0 erupt_start],[4 4],'LineWidth',3,'Color',cmap(2,:),'LineStyle','--');
plot([erupt_end max(time)],[6.1 6.1],'LineWidth',3,'Color',cmap(3,:),'LineStyle','--');
h = plot(time(3:end), Q_filt(3:end),'k'); 
h.LineWidth = 2;
ylabel('Quality Factor');
xlim([min(t_hrs) max(t_hrs)]);
h.Parent.XTick = [0 12 24 36 48 60 72];
h.Parent.XTickLabel = {};
ylim([0 18]);
box on;

% plot eruption time window
h = rectangle('Position',[erupt_start 0 3 18]);
h.EdgeColor = [1 0 0 0.5];
h.FaceColor = [1 0 0 0.5];


%% spectra before and after %%

% spectra before onset of flank eruption
subplot(4,2,7);
Data_before_stack = mean(abs(Data_before_norm_save),2);
[f0_before, Q_before] = resPeakProps(f,abs(Data_before_stack(1:NFFT/2+1)));
f0_before = 0.60; % average from hourly values
plot([f0_before f0_before],[0 1],'LineWidth',3,'Color',cmap(2,:),'LineStyle','--')
hold on;
plot(f, abs(Data_before_stack(1:NFFT/2+1))./...
    max(abs(Data_before_stack(1:NFFT/2+1)))','Color','k','LineWidth',2);
xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
xlim([0 3])
hold on;


% spectra after onset of flank eruption
subplot(4,2,8);
Data_after_stack = mean(abs(Data_after_norm_save),2);
[f0_after, Q_after] = resPeakProps(f,abs(Data_after_stack(1:NFFT/2+1)));
f0_after = 0.36; % average from hourly values 
plot([f0_after f0_after],[0 1],'LineWidth',3,'Color',cmap(3,:),'LineStyle','--')
hold on;
plot(f, abs(Data_after_stack(1:NFFT/2+1))./...
    max(abs(Data_after_stack(1:NFFT/2+1)))','Color','k','LineWidth',2);
xlabel('Frequency (Hz)');
%ylabel('Normalized Amplitude');
xlim([0 3])
hold on;



