%% MAKE FIG COMBINED %%
%
% Make combined figure showing both of the numerical verification exercises
% (for narrow and wide conical horn)

clear all; clc; %close all;
path(pathdef); % clear path

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',16);
cmap = get(gca,'ColorOrder');
alpha = 0.5;

% add paths
addpath ../source/SBPoperators/
addpath ../source/resonance/
addpath Functions/
addpath Data/Fig5_NumericalBenchmarking/

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 1100 600]);

%% NARROW CRATER %%

%%%%%%%%%%%%%%%%%%
%%% INFRA FDTD %%%
%%%%%%%%%%%%%%%%%%


%%% excess pressure in time domain %%%
sta_str = 'STA4_pressure.txt';
sta = load(sta_str);
t = sta(:,1); % station time
p = sta(:,2); % station excess pressure
amp = 1000; % amplitude factor to scale pressure by
p = p*amp;

% plot excess pressure in time domain
subplot(2,6,[3 4]);
h1 = plot(t, p,'Color',cmap(1,:));
h1.Color(4) = alpha;
hold on;
xlabel('Time (s)');
ylabel('\Delta p');

%%% excess pressure in frequency domain %%%
dt = t(2)-t(1); % time step
P = fft(p)*dt; % fourier transform
N = length(t); % number of time samples
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
f = omega/(2*pi); % frequency vector

% plot excess pressure in frequecy domain
subplot(2,6,[5 6]);
h1 = plot(f(1:floor(N/2)), abs(P(1:floor(N/2))));
h1.Color(4) = alpha;
xlim([0 2]);
set(h1.Parent,'YTick',[]);
hold on;
xlabel('Freq (Hz)');

%%%%%%%%%%%%%%
%%% RES 1D %%%
%%%%%%%%%%%%%%

%%% set parameters %%%
gamma = 1; % isothermal to match analytical solutions
r = 200; % distance from crater to station
M = problemParameters_validation(gamma, r); % load problem parameters
M.rhoA = 1.22;
M.cA = 340;
M.rhoC = 1.22;
M.cC = 340;
M.pC = M.cC^2*M.rhoC/M.gamma;
style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method
%%% geometry %%%
shape = load('res1d_geometry.txt');
depth = shape(1,1);

%%% plot geometry %%%
subplot(2,6,[1 2]);
h1 = plot(shape(:,2), shape(:,1),'Color',cmap(1,:)); h1.Color(4) = alpha;
set(gca,'YDir','Reverse');
hold on;
h1 = plot(-shape(:,2), shape(:,1),'Color',cmap(1,:)); h1.Color(4) = alpha;
plot([-shape(1,2) shape(1,2)],[shape(1,1) shape(1,1)],'r');
h1 = plot([shape(end,2) 450],[0 0],'Color',cmap(1,:)); h1.Color(4) = alpha;
h1 = plot([-shape(end,2) -350],[0 0],'Color',cmap(1,:)); h1.Color(4) = alpha;
xlim([-250 450]);
ylim([-50 200]);
xlabel('Radius (m)');
ylabel('Depth (m)');
h1 = plot(200,-10,'v','MarkerFaceColor',cmap(1,:),'MarkerEdgeColor',cmap(1,:),'MarkerSize',12);
h1.Color(4) = alpha;

%%% load source %%%
src = load('monopole_src_1.txt');
ts = src(:,1)';
s = src(:,2)';
s = s*amp;

% source from infraFDTD is in terms of a mass flux. Source for res1D is in
% terms of a volume flux. Therefore, divide by density of air and radius at
% base of crater
s = s./M.rhoA;
area = pi*shape(1,2)^2;

ts = ts(1:end-1); % formulas assume that there are an even number of time samples
s = s(1:end-1);

T = max(ts);
N = length(ts);
dt = T/N;

Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples


S = fft(s/area)*dt;
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
f = omega/(2*pi); % frequency vector

% plot source in time domain
subplot(2,6,7);% figure(1); 
plot(ts,s,'k');
title('Source');
xlabel('Time (s)');
ylabel('m^3/s');
xlim([0 5]);

% plot source in frequency domain
subplot(2,6,8);
h = plot(f(1:N/2+1), abs(S(1:N/2+1))./max(abs(S(1:N/2+1))),'k');
set(h.Parent,'YTick',[]);
xlim([0 2]);
xlabel('Freq (Hz)');

%%% COMPUTE GREEN'S FUNCTION %%%
tic
A = resonance1d(shape, depth, freq, Nf, style, order, M);

% convolve Green's function with source
pres = A.P(1:N/2+1) .* S(1:N/2+1);
toc



% plot infrasound signal in frequency domain
subplot(2,6,[5 6]);
plot(A.f(1:N/2+1), abs(pres),'Color',cmap(1,:),'LineStyle',':');

% invert infrasound signal to time domain
presF = [pres conj(pres(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
presF(N/2+1) = real(presF(N/2+1)); % entry at Nyquist must be real
presTime = ifft(presF,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain

% plot infrasound signal in time domain
subplot(2,6,[3 4]);
plot(ts, presTime,'Color',cmap(1,:),'LineStyle',':');


%% WIDE CRATER %%

clear all
cmap = get(gca,'ColorOrder');
alpha = 0.5;

%%%%%%%%%%%%%%%%%%
%%% INFRA FDTD %%%
%%%%%%%%%%%%%%%%%%

%%% excess pressure in time domain %%%
sta_str = 'STA4_pressure_wide.txt';
sta = load(sta_str);
t = sta(:,1); % station time
p = sta(:,2); % station excess pressure
amp = 1000; % amplitude factor to scale pressure by
p = p*amp;

% plot excess pressure in time domain
subplot(2,6,[9 10]);
h1 = plot(t, p,'Color',cmap(2,:));
h1.Color(4) = alpha;
hold on;
xlabel('Time (s)');
ylabel('\Delta p');
xlim([0 10])

%%% excess pressure in frequency domain %%%
dt = t(2)-t(1); % time step
P = fft(p)*dt; % fourier transform
N = length(t); % number of time samples
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
f = omega/(2*pi); % frequency vector

% plot excess pressure in frequecy domain
subplot(2,6,[11 12]);
h1 = plot(f(1:floor(N/2)), abs(P(1:floor(N/2))),'Color',cmap(2,:));
h1.Color(4) = alpha;
xlim([0 2]);
set(h1.Parent,'YTick',[]);
hold on;
xlabel('Freq (Hz)');

%%%%%%%%%%%%%%
%%% RES 1D %%%
%%%%%%%%%%%%%%

%%% set parameters %%%
gamma = 1; % isothermal to match analytical solutions
r = 400; % distance from crater to station
M = problemParameters_validation(gamma, r); % load problem parameters
M.rhoA = 1.22;
M.cA = 340;
M.rhoC = 1.22;
M.cC = 340;
M.pC = M.cC^2*M.rhoC/M.gamma;
style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method

%%% geometry %%%
shape = load('res1d_geometry_wide.txt');
depth = shape(1,1);

%%% geometry %%%
subplot(2,6,[1 2]);
h1 = plot(shape(:,2), shape(:,1),'Color',cmap(2,:)); h1.Color(4) = alpha;
set(gca,'YDir','Reverse');
hold on;
h1 = plot(-shape(:,2), shape(:,1),'Color',cmap(2,:)); h1.Color(4) = alpha;
plot([-shape(1,2) shape(1,2)],[shape(1,1) shape(1,1)],'r');
h1 = plot([shape(end,2) 450],[0 0],'Color',cmap(2,:)); h1.Color(4) = alpha;
h1 = plot([-shape(end,2) -350],[0 0],'Color',cmap(2,:)); h1.Color(4) = alpha;
xlim([-350 450]);
ylim([-50 200]);
xlabel('Radius (m)');
ylabel('Depth (m)');
h1 = plot(400,-10,'v','MarkerFaceColor',cmap(2,:),'MarkerEdgeColor',cmap(2,:),...
    'MarkerSize',12); h1.Color(4) = alpha;


%%% load source %%%
src = load('monopole_src_1.txt');
ts = src(:,1)';
s = src(:,2)';
s = s*amp;

% source from infraFDTD is in terms of a mass flux. Source for res1D is in
% terms of a volume flux. Therefore, divide by density of air and radius at
% base of crater
s = s./M.rhoA;
area = pi*shape(1,2)^2;

ts = ts(1:end-1); % formulas assume that there are an even number of time samples
s = s(1:end-1);

T = max(ts);
N = length(ts);
dt = T/N;

Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples


S = fft(s/area)*dt;
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
f = omega/(2*pi); % frequency vector

% plot source in time domain
subplot(2,6,7);% figure(1); 
plot(ts,s,'k');
title('Source');
xlabel('Time (s)');
ylabel('m^3/s');
xlim([0 5]);

% plot source in frequency domain
subplot(2,6,8);
h = plot(f(1:N/2+1), abs(S(1:N/2+1))./max(abs(S(1:N/2+1))),'k');
set(h.Parent,'YTick',[]);
xlim([0 2]);
xlabel('Freq (Hz)');

%%% COMPUTE GREEN'S FUNCTION %%%
A = resonance1d(shape, depth, freq, Nf, style, order, M);

% convolve Green's function with source
pres = A.P(1:N/2+1) .* S(1:N/2+1);

% plot infrasound signal in frequency domain
subplot(2,6,[11 12]);
plot(A.f(1:N/2+1), abs(pres),'Color',cmap(2,:),'LineStyle',':');

% invert infrasound signal to time domain
presF = [pres conj(pres(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
presF(N/2+1) = real(presF(N/2+1)); % entry at Nyquist must be real
presTime = ifft(presF,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain

% plot infrasound signal in time domain
subplot(2,6,[9 10]);
plot(ts, presTime,'Color',cmap(2,:),'LineStyle',':');