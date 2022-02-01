% FORWARD %

clear all; clc;
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',18);
cmap = get(gca,'ColorOrder');

addpath ../source/inv/
addpath ../source/resonance/
addpath ../source/SBPoperators/

%% TOPOGRAPHY %%
%
% Plot topographic profiles from DEM's provided by INGV
% These were emailed by Emanuela De Beni on October 29 and 30.
% Cross-sections were taken from the 2021 March 3 DEM

W = readmatrix('w_section.txt');

% shift profile so that outlet height is zero and crater is centered at zero
z_shift = max(W(:,2));
r_shift = 162;

W0(:,1) = W(:,1)-r_shift;
W0(:,2) = -(W(:,2)-z_shift);

%% INITIALIZE CRES %%

% set the parameters for the resonance1d calculations
sps = 20; % samples per second
dt = 1/sps; % time step
N = 1000; % number of grid points (formulas assume even N)
T = dt*N; % total time (s)
Nyquist = 1/(2*dt); % Nyquist frequency (Hz)
Nf = N/2+1; % number of frequency samples
freq = [0 Nyquist]; % frequency range (Hz)
discrParams = [T N Nf Nyquist dt]; % save parameters into array

order = 4; % order of numerical scheme (4, 6 or 8)
style = 'baffled piston'; % acoustic radiation model ('monopole' or ' baffled piston')

% filtering properties
filterband = [0.5 Nyquist-0.1];
filterorder = 4; % order of butterworth filter
Fs = sps; % sampling frequency (same as samples per second)
filterProps = [filterband, filterorder, Fs]; % filter properties - same as for data (except data is high pass but this is a band filter)

% geometry
load EtnaGeomCRes.mat

figure(1); clf;
plot(shape(:,2), shape(:,1),'Color',cmap(1,:));
hold on;
set(gca,'YDir','Reverse');
plot(-shape(:,2), shape(:,1),'Color',cmap(1,:));
ylim([-50 300]); grid on; hold on;
xlabel('Radius (m)');
ylabel('Depth (m)');

% plot DEM profile
W = readmatrix('w_section.txt');

% shift profile so that outlet height is zero and crater is centered at zero
z_shift = max(W(:,2))-10;
r_shift = 160;

W0(:,1) = W(:,1)-r_shift;
W0(:,2) = -(W(:,2)-z_shift);

plot(W0(:,1), W0(:,2),'k');

% source
srcParams = 0.3;
srcStyle = 'Brune';
discrParamsSrc = [T N];
[S, f, s, ts] = sourceFunction(1, srcParams, srcStyle, discrParamsSrc);

figure(2); clf;
subplot(2,1,1);
plot(ts, s./max(s)); xlim([10 20]); hold on;
xlabel('Time (s)');
ylabel('Normalized Amplitude');
title('Source: Time Domain');
subplot(2,1,2);
plot(f(1:N/2+1), abs(S(1:N/2+1))./max(abs(S(1:N/2+1)))); hold on;
xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
title('Source: Time Domain');

% temperature
depth = 300;
T1 = 1100;
T2 = 1000;
T3 = 450;
D1 = 55;
Ta = -7.5;

idx_z_crater = D1+1;
idx_z_conduit = depth-D1;

T_profile1 = linspace(T3,T2,idx_z_crater);
T_profile2 = linspace(T2, T1, idx_z_conduit);
T_profile = [T_profile1 T_profile2]';
Tp = smoothdata(T_profile,'lowess');
Z = [0:1:depth]';
Z = flipud(Z);
Tp = flipud(Tp);
figure(3); clf;
plot(Tp,Z); hold on;
set(gca,'YDir','Reverse');
ylabel('Depth (m)');
xlabel('Temperature (C)');
grid on;
M = problemParameters_temp(Tp,Ta,Z);
   

%% CRES  %%

% CRes outputs
DEPTH_VECTOR = depth:-10:50;
L = length(DEPTH_VECTOR);
SPECTRA = zeros(N/2+1, L,1);
SPECTRA_NORM = zeros(N/2+1, L,1);

% CRes inversion

%%% CRes %%%
tic
f0_sim = zeros(L,5);
for j = 1:L
    DEPTH_VECTOR(j)
    
    res = resonance1d_temp(shape, DEPTH_VECTOR(j), freq, Nf, style, order, M); % compute transfer function
    
    sim.f = res.f; % frequency vector
    sim.P = (res.P(1:N/2+1).*S(1:N/2+1));
    sim.PNorm = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1));
    
    SPECTRA(:,j) = abs(sim.P);
    SPECTRA_NORM(:,j) = abs(sim.PNorm);
    %[f0_sim(j), Q_sim(j)] = resPeakProps(sim.f,abs(sim.P));
    [xmax,imax,xmin,imin] = extrema(abs(sim.P));
    
    for i = 1:min([5,length(imax)])
        f0_sim(j,i) = sim.f(imax(i));
        h = vline(sim.f(imax(i)));
    end
end
f = sim.f;
toc

%% save outputs

save('CRES.mat','f0_sim','DEPTH_VECTOR','SPECTRA','SPECTRA_NORM','f');

%% FORMAT PEAK FREQ VECTOR %%

figure(4); clf;
plot(f0_sim(:,1), DEPTH_VECTOR);
set(gca,'YDir','Reverse');
hold on; grid on;
xlabel('Peak Freq. (Hz)');
ylabel('Depth (m)');
plot(f0_sim(:,2), DEPTH_VECTOR);

idx1 = find(DEPTH_VECTOR == 106);
idx2 = find(DEPTH_VECTOR == 105);

plot(f0_sim(1:idx1,1), DEPTH_VECTOR(1:idx1),'Color',cmap(3,:));
plot(f0_sim(idx2:end,2), DEPTH_VECTOR(idx2:end),'Color',cmap(3,:));
plot(f0_sim(idx2:end,3), DEPTH_VECTOR(idx2:end),'Color',cmap(4,:));

idx3 = find(DEPTH_VECTOR == 86);
idx4 = find(DEPTH_VECTOR == 66);

F0 = [f0_sim(1:idx1,1); f0_sim(idx2:end,2)];
F0(idx3) = f0_sim(idx3,1);
F0(idx4) = f0_sim(idx4,3);
D = DEPTH_VECTOR;
plot(F0, D,'k');
hold on;
F0_raw = F0;
F0 = smoothdata(F0,'sgolay');
plot(F0,D);

%% CRes outputs %%%
figure(5); clf; 
plot(F0, D);
set(gca,'YDir','Reverse');
hold on; grid on;
xlabel('Peak Freq. (Hz)');
ylabel('Depth (m)');

% inversion values
peak_freq = F0;
d_sim = D; 

% spectra
[FF,DD] = meshgrid(f,DEPTH_VECTOR);
figure(6); clf; subplot(1,2,1);
surf(FF,DD,SPECTRA_NORM');
view(2);
shading interp
ylabel('Depth (m)');
xlabel('Frequency (Hz)');
colorbar
xlim([0 5]);
ylim([min(DEPTH_VECTOR) max(DEPTH_VECTOR)])
set(gca,'YDir','Reverse');
title('Normalized Spectra');
hold on;
plot3(F0,D,10*ones(size(D)),'r');

subplot(1,2,2);
surf(FF,DD,SPECTRA');
view(2);
shading interp
ylabel('Depth (m)');
xlabel('Frequency (Hz)');
colorbar
xlim([0 5]);
ylim([min(DEPTH_VECTOR) max(DEPTH_VECTOR)])
set(gca,'YDir','Reverse');
title('Spectra');

%% save outputs %%

% save('forwardOutput.mat');







