%% Inversion Ascent%%
% Calculate magma ascent velocity/volumetric rate for given inverted depth

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
cmap = get(gca,'ColorOrder');

%% load depth
load inversionDepth_200C.mat

t = out1.t-24;
D0 = out1.D;

figure(1); clf;
plot(t, D0,'.');
xlabel('Time (hours0');
ylabel('Depth (m)');
set(gca,'YDir','Reverse');

% remove depths less than 200 m
D = D0;
for j = 1:5
    for i = 1:length(t)
        if D(j,i) > 190
            D(j,i) = NaN;
        end
    end
end

figure(2); clf;
plot(t, D,'.');




%% smooth depth
% divide into two time sections
tidx = find(t > 2,1,'first');
tidx2 = find(t > 19,1,'first');
t1 = t(1:tidx-1);
t2 = t(tidx:tidx2-1);
t3 = t(tidx2:end);


D1 = D(:,1:tidx-1);
D2 = D(:,tidx:tidx2-1);
D3 = D(:,tidx2:end);

sl_time = 2; % smoothing length
dt = t(2)-t(1);
sl = round(sl_time/dt);

figure(2); clf;
subplot(3,1,1);
figure(3); clf;
for j = 1:5
    DS1(j,:) = movmean(D1(j,:),sl,'omitnan');
    DS2(j,:) = movmean(D2(j,:),sl,'omitnan');
    DS3(j,:) = movmean(D3(j,:),sl,'omitnan');
    figure(3);
    subplot(3,2,j);
    plot(t1,D1(j,:),'.','Color',cmap(j,:));
    hold on;
    plot(t2,D2(j,:),'.','Color',cmap(j,:));
    plot(t3,D3(j,:),'.','Color',cmap(j,:));

    plot(t1,DS1(j,:),'-','Color',cmap(j+1,:));
    hold on;
    plot(t2,DS2(j,:),'-','Color',cmap(j+1,:));
    plot(t3,DS3(j,:),'-','Color',cmap(j+1,:));
end
xlabel('Time (hours0');
ylabel('Depth (m)');
set(gca,'YDir','Reverse');



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
R3 = zeros(size(DS3));
A1 = zeros(size(DS1));
A2 = zeros(size(DS2));
A3 = zeros(size(DS3));
for i = 1:5
    
    % plot depth
    subplot(1,3,1);
    plot(t1, DS1(i,:),'.','Color',cmap(i,:));
    hold on;
    plot(t2, DS2(i,:),'.','Color',cmap(i,:));
    plot(t3, DS3(i,:),'.','Color',cmap(i,:));
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
    
    for j = 1:length(DS3(i,:))
        d_find = DS3(i,j);
        if isnan(d_find)
            R3(i,j) = NaN;
        else
            d_idx = find(d_find >= D,1,'first');
            r_tmp = R0(d_idx);
            R3(i,j) = r_tmp;
        end
    end
    
    % plot radius
    subplot(1,3,2);
    plot(t1, R1(i,:),'.','Color',cmap(i,:));
    hold on;
    plot(t2, R2(i,:),'.','Color',cmap(i,:));
    plot(t3, R3(i,:),'.','Color',cmap(i,:));
    title('Radius');
    
    % plot cross-sectional area
    A1(i,:) = pi*R1(i,:).^2;
    A2(i,:) = pi*R2(i,:).^2;
    A3(i,:) = pi*R3(i,:).^2;
    subplot(1,3,3);
    plot(t1, A1(i,:),'.','Color',cmap(i,:));
    hold on;
    plot(t2, A2(i,:),'.','Color',cmap(i,:));
    plot(t3, A3(i,:),'.','Color',cmap(i,:));
    title('Area');
end


%% Differentiate to get magma ascent rate %%

t1d = t1(1:end-1)+dt/2; % time vectors for differentiation
t2d = t2(1:end-1)+dt/2; % time vectors for differentiation
t3d = t3(1:end-1)+dt/2; % time vectors for differentiation

figure(10); clf; subplot(2,1,1);

for i = 1:5
    M1(i,:) = diff(DS1(i,:))./dt; % magma ascent rate in m/hr
    M2(i,:) = diff(DS2(i,:))./dt; % magma ascent rate in m/hr
    M3(i,:) = diff(DS3(i,:))./dt; % magma ascent rate in m/hr
    
    plot(t1d, -M1(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    hold on;
    plot(t2d, -M2(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    plot(t3d, -M3(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    ylim([-10 50]);
    grid on;

end
xlabel('Time (hours)');
ylabel('Ascent Velocity (m/hr)');
%ylim([-10 50]);
grid on;


%% Calculate volume flux of magma

% calculate volume flux by multiplying cross-sectional area with magma
% ascent velocity
A1d = A1(:,1:end-1);
A2d = A2(:,1:end-1);
A3d = A3(:,1:end-1);
VolFlux1 = zeros(size(A1d));
VolFlux2 = zeros(size(A2d));
VolFlux3 = zeros(size(A3d));
%figure(8); clf;
figure(10); subplot(2,1,2);
for i = 1:5
    VolFlux1(i,:) = -M1(i,:) .* A1d(i,:);
    VolFlux2(i,:) = -M2(i,:) .* A2d(i,:);
    VolFlux3(i,:) = -M3(i,:) .* A3d(i,:);
    plot(t1d, VolFlux1(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    hold on;
    plot(t2d, VolFlux2(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    plot(t3d, VolFlux3(i,:),'.','Color',cmap(i,:),'MarkerSize',10);
    ylim([-5000 15000])
grid on;
end
grid on;
xlabel('Time (hr since 00:00 19 Feb 2021)');
ylabel('Magma Volume Flux (m^3/hr)');

time = [t1d t2d t3d];
vel = [-M1 -M2 -M3];
vol = [VolFlux1 VolFlux2 VolFlux3];

save('inversionAscent.mat','time','vel','vol');


