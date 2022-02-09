% BACKWARD %

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
cmap = get(gca,'ColorOrder');

addpath ../source/inv/
addpath ../source/resonance/
addpath ../source/SBPoperators/
 
%% LOAD INFRASOUND DATA %%

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



%% Inversion %%%

load forwardOutput.mat

DEPTH_F0 = zeros(5,N);

% iterate over stations
figure(2); clf;
for j = 1:5
    j
    
    % iterate over time
    for k = 1:length(t)
        
        fval = f0smooth(j,k); % peak frequency value from data at station j at time k
        if isnan(fval)
            DEPTH_F0(j,k) = NaN;
        else
            tmp_f0 = abs(peak_freq - fval); % subtract freq value of interest from simulations
            [~,~,tmpmin_f0,imin_f0] = extrema(tmp_f0); % find local minimun
            
            if length(imin_f0) == 1
                DEPTH_F0(j,k) = d_sim(imin_f0);
            else
                DEPTH_F0(j,k) = d_sim(imin_f0(1));
            end
        end
        
        if DEPTH_F0(j,k) == depth
            DEPTH_F0(j,k) = NaN;
        end
        
    end
    
    % plot results. each station should be a different plot with the
    % three different sound speeds
    plot(t, DEPTH_F0(j,:),'.','MarkerSize',10);
    hold on;
    xlabel('Time (hrs) since 00:00 on Feb 19, 2021');
    ylabel('Depth (m)');
    set(gca,'YDir','Reverse');
    grid on;
    xlim([min(tplot) max(tplot)]);
    ylim([50 300])
    title('Non-smoothed Depths');
    
end

return


%% Divide data into the two time sections. Smooth data using movmean

tidx = find(t > 30,1,'first');

t1 = t(1:tidx-1);
t2 = t(tidx:end);

D1 = DEPTH_F0(:,1:tidx-1);
D2 = DEPTH_F0(:,tidx:end);

sl_time = 0.75; % one hour smoothing length
dt = t(2)-t(1);
sl = round(sl_time/dt)

figure(3); clf;
for j = 1:5
    DS1(j,:) = movmean(D1(j,:),sl);
    DS2(j,:) = movmean(D2(j,:),sl);
    
    plot(t1,D1(j,:),'.','Color',cmap(j,:));
    hold on;
    plot(t2,D2(j,:),'.','Color',cmap(j,:));


    plot(t1,DS1(j,:),'-','Color',cmap(j,:));
    hold on;
    plot(t2,DS2(j,:),'-','Color',cmap(j,:));
end

set(gca,'YDir','Reverse');
xlabel('Time');
ylabel('Depth');
title('Smoothed Depths');

%% Calculate radius and area as a function of time %%

% load and plot crater geometry
load EtnaGeomCRes.mat
R0 = shape(:,2);
D = shape(:,1);
figure(4); clf;
plot(R0,D,'k');
hold on;
set(gca,'YDir','Reverse');
plot(-R0,D,'k');
xlabel('Radius (m)');
ylabel('Depth (m)');

% calculate cross-sectional area as a function of depth
A = pi*R0.^2;
figure(5); clf;
plot(A,D);
xlabel('Cross-Sectional Area (m^2)');
ylabel('Depth (m)');
set(gca,'YDir','Reverse');


% plot depth, radius, and cross-sectional area as a function of time
figure(6); clf;
R1 = zeros(size(DS1));
R2 = zeros(size(DS2));
A1 = zeros(size(DS1));
A2 = zeros(size(DS2));
for i = 1:5
    
    % plot depth
    subplot(1,3,1);
    plot(t1, DS1(i,:),'.','Color',cmap(i,:));
    hold on;
    plot(t2, DS2(i,:),'.','Color',cmap(i,:));
    set(gca,'YDir','Reverse');
    title('Depth');

    % for each depth, find corresponding value of radii
    for j = 1:length(DS1(i,:))
        d_find = DS1(i,j);
        if isnan(d_find)
            R1(i,j) = NaN;
        else
            d_idx = find(D <= d_find,1,'first')
            r_tmp = R0(d_idx);
            R1(i,j) = r_tmp;
        end
    end
    
    for j = 1:length(DS2(i,:))
        d_find = DS2(i,j);
        if isnan(d_find)
            R2(i,j) = NaN;
        else
            d_idx = find(d_find >= D,1,'first');
            r_tmp = R0(d_idx);
            R2(i,j) = r_tmp;
        end
    end
    
    % plot radius
    subplot(1,3,2);
    plot(t1, R1(i,:),'.','Color',cmap(i,:));
    hold on;
    plot(t2, R2(i,:),'.','Color',cmap(i,:));
    title('Radius');
    
    % plot cross-sectional area
    A1(i,:) = pi*R1(i,:).^2;
    A2(i,:) = pi*R2(i,:).^2;
    subplot(1,3,3);
    plot(t1, A1(i,:),'.','Color',cmap(i,:));
    hold on;
    plot(t2, A2(i,:),'.','Color',cmap(i,:));
    title('Area');
end


%% Differentiate to get magma ascent rate %%

t1d = t1(1:end-1)+dt/2; % time vectors for differentiation
t2d = t2(1:end-1)+dt/2; % time vectors for differentiation

figure(7); clf;

for i = 1:5
    M1(i,:) = diff(DS1(i,:))./dt; % magma ascent rate in m/hr
    M2(i,:) = diff(DS2(i,:))./dt; % magma ascent rate in m/hr
    
    plot(t1d, M1(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    hold on;
    plot(t2d, M2(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    ylim([-10 50]);
    grid on;

end
xlabel('Time (hours)');
ylabel('Ascent Velocity (m/hr)');
ylim([-10 50]);
grid on;


%% Calculate volume flux of magma

% calculate volume flux by multiplying cross-sectional area with magma
% ascent velocity
A1d = A1(:,1:end-1);
A2d = A2(:,1:end-1);
VolFlux1 = zeros(size(A1d));
VolFlux2 = zeros(size(A2d));
figure(8); clf;
for i = 1:5
    VolFlux1(i,:) = -M1(i,:) .* A1d(i,:);
    VolFlux2(i,:) = -M2(i,:) .* A2d(i,:);
    plot(t1d, VolFlux1(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    hold on;
    plot(t2d, VolFlux2(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    ylim([-5000 15000])
grid on;
end
grid on;
xlabel('Time (hr since 00:00 19 Feb 2021)');
ylabel('Magma Volume Flux (m^3/hr)');

%% save %%

save('inversionOutput.mat');

%% Make figure for paper

figure(10); clf;
for i = 1:5
    
    % frequency
    subplot(4,1,1);
    plot(t,f0smooth(i,:),'.','MarkerSize',10);
    ylabel('Frequency (Hz)');
    grid on; hold on;
      
    % depth
    subplot(4,1,2);
    plot(t1,D1(i,:),'.','Color',cmap(i,:));
    hold on;
    plot(t2,D2(i,:),'.','Color',cmap(i,:));
    set(gca,'YDir','Reverse');
    ylabel('Depth (m)'); grid on;
    
    % ascent velocity
    subplot(4,1,3);
    plot(t1d, -M1(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    hold on; grid on;
    plot(t2d, -M2(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    ylabel('Ascent Velocity (m/hr)');
    ylim([-20 60])
    
    % volume flux
    subplot(4,1,4);
    plot(t1d, VolFlux1(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    hold on; grid on;
    plot(t2d, VolFlux2(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    ylabel('Volume Flux');
    ylim([-5000 15000])
end





