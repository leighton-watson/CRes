%% format data

clear all; clc;

load etna2021_data_f0_Q.mat;

% isolate data of interest when signal from volcano dominates
tplot = [20 24 30 48];

% time vector in hours
t = data_t/3600;
f0 = data_f0;

% remove entries outside time windows of interest
tidx1 = find(t>tplot(1),1,'first');
t(1:tidx1) = [];
f0(:,1:tidx1) = [];

tidx2 = find(t>tplot(2),1,'first');
tidx3 = find(t>tplot(3),1,'first');
t(tidx2:tidx3) = [];
f0(:,tidx2:tidx3) = [];

tidx4 = find(t>tplot(4),1,'first');
t(tidx4:end) = [];
f0(:,tidx4:end) = [];

% apply smoothing
f0smooth = zeros(size(f0));
filt_length = 10;
for i = 1:5
    f0smooth(i,:) = medfilt1(f0(i,:),filt_length);
end

% plot peak frequency
figure(1); clf;
for j = 1:5
    plot(t,f0smooth(j,:),'.','MarkerSize',10);
    ylabel('Frequency (Hz)');
    xlim([min(tplot) max(tplot)]);
    grid on;
    hold on;
      
    F_MAX(j) = max(f0smooth(j,:));
end

FMAX = max(F_MAX);

save('etna2021_data_peak_freq.mat','t','f0smooth');