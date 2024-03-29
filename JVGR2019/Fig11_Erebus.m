%% MAKE FIG 2 %%
%
% Make figure comparing infraFDTD and res1D for Erebus
% Produce plots as separate figures that can be edited in Illustrator


clear all; clc; %close all;
path(pathdef); % clear path

set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',20);
cmap = get(gca,'ColorOrder');

% add paths
addpath ../source/SBPoperators/
addpath ../source/resonance/
addpath Data/Fig11_Erebus/


%% CONTOUR PLOT %%

elev = load('sur_coords.txt');

x = elev(:,1); 
y = elev(:,2); 
z = elev(:,3);

nx = 499; ny = 499; nz = 269;

X = reshape(x,nx,ny);
Y = reshape(y,nx,ny);
Z = reshape(z,nx,ny);

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 900 600 500]);

h1 = contour(X,Y,Z,25);
shading interp
axis tight
%demcmap([2740 2850]);
demcmap([min(min(Z)) max(max(Z))]);
hold on;
xlabel('X (m)');
ylabel('Y (m)');
h2 = colorbar;
ylabel(h2,'Elevation (m a.s.l.)');
%set(h2,'Ticks',[2750 2775 2800 2825 2850]);
view(2);

% station position
num_sta = 6;
y_sta = [1550 1600 1780 1940 1025 400];
x_sta = [1085 1260 1210 1565 1140 1075];
z_sta = [3414 3472 3421 3404 3439 3271]+10;

for j = 1:num_sta
    plot3(x_sta(j), y_sta(j), z_sta(j), 'v','MarkerSize',12,...
        'MarkerEdgeColor',cmap(j,:),'MarkerFaceColor',cmap(j,:));
end

% source position
y_src = 1650;
x_src = 1145;
z_src = 3420;

plot3(x_src, y_src, z_src, 'o','MarkerSize',12,...
    'MarkerEdgeColor','k','MarkerFaceColor','k');



%% CRATER PROFILE %%
 
figHand2 = figure(2); clf;
set(figHand2,'Position',[100 300 600 500]);

load ErebusDEM_Cyl.mat;   % load DEM

% extract cross-section in radial slices through DEM
theta_multi = linspace(0,2*pi,100);
r = 25:150;
rz_save = zeros(length(theta_multi), length(r));
alpha = 0.7;

% plot radial slices through 3D geometry
for i = 1:length(theta_multi)
    i
    theta = theta_multi(i);
    rx = sin(theta).*r;
    ry = cos(theta).*r;

    rz = interp2(X,Y,Z,rx,ry); % interpolate onto regular radius vector
    plot(r,rz,'color',[0 0 0]+alpha,'LineWidth',10); hold on;
    plot(-r,rz,'color',[0 0 0]+alpha,'LineWidth',10);
    
    rz_save(i,:) = rz; % save results
end
rz_avg = mean(rz_save); % compute mean radius for axisymmetric crater
 
 
 
% axisymmetric res1D geometry
h3 = plot(r, rz_avg,'k');
plot(-r, rz_avg,'k');

% plot lava lake
plot([-25 25],[min(rz_avg) min(rz_avg)],'r')

% format figure
xlim([-150 150])
ylim([min(rz_avg) 3750]);
ylabel('Elevation (m a.s.l.)');
xlabel('Radius (m)');
%set(h3.Parent,'YTick',[2750 2775 2800 2825 2850]);



%% SOURCE %%

figHand3 = figure(3); clf;
set(figHand3,'Position',[100 300 600 300]);

%%% source in time domain %%%
src = load('monopole_src_1.txt');
ts = src(:,1); % source time
s = src(:,2); % source amplitude 
amp = 1000; % factor to scale amplitude of source by
s = s*amp; % increase amplitude of source

subplot(1,2,1);
%figHand3 = figure(3); clf;
%set(figHand3,'Position',[100 100 600 250]);
%subplot(1,2,2);
plot(ts, s./1.22,'k');
xlabel('Time (s)');
ylabel('m^3/s');
hold on;
xlim([0 3])

%%% source in frequency domain %%%
dt = ts(2)-ts(1); % time step
S = fft(s)*dt; % fourier transform 
N = length(ts); % number of time samples
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
f = omega/(2*pi); % frequency vector

subplot(1,2,2);
%h = plot(f(1:floor(N/2)), abs(S(1:floor(N/2)))./max(abs(S(1:floor(N/2)))));
h = plot(f(1:floor(N/2)), abs(S(1:floor(N/2))),'k');
%set(h.Parent,'YTick',[]);
xlim([0 3]);
ylim([0 650])
xlabel('Freq (Hz)');
ylabel('(m^3/s)/Hz');
hold on;

%% INFRA FDTD %%

figHand4 = figure(4); clf;
set(figHand4,'Position',[600 100 500 1500]);

figHand5 = figure(5); clf;
set(figHand5,'Position',[600 100 500 1500]);

sta_str1 = 'ST';
sta_str2 = '_pressure.txt';
num_sta = 6;
for i = 1:num_sta
    
    % load data
    sta_str = strcat(sta_str1,num2str(i),sta_str2);
    sta = load(sta_str);
    
    %%% excess pressure in time domain
    t = sta(:,1); % station time
    p = sta(:,2); % station excess pressure
    p = p*amp; % scale amplitude by the same amount as the source was scaled
    
    figure(4);
    subplot(num_sta,1,i);
    h1 = plot(t,p,'Color',cmap(i,:));
    h1.Color(4) = 0.7;
    hold on;
    %xlim([0 10])
    
    p_3D(i) = max(p);

    %%% excess pressure in frequency domain %%%
    dt = t(2)-t(1); % time step
    P = fft(p)*dt; % fourier transform
    N = length(t); % number of time samples
    omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
    f = omega/(2*pi); % frequency vector
    
    figure(5);
    subplot(num_sta,1,i);
    h2 = plot(f(1:floor(N/2)), abs(P(1:floor(N/2))),'Color',cmap(i,:));
    h2.Color(4) = 0.7;
    hold on;
    xlim([0 3])
    
    
end 


%% RES 1D


rsta_multi = sqrt((x_sta-x_src).^2 + (y_sta-y_src).^2);

%%% set parameters %%%
gamma = 1; % isothermal to match infraFDTD solutions
style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method

M.gamma = gamma;
M.R = 287.06;
M.rhoA = 1.22;
M.cA = 340;
M.rhoC = 1.22;
M.cC = 340;
M.pC = M.cC^2*M.rhoC/M.gamma;
M.TA = M.cA.^2/(M.gamma*M.R) - 273.15;
M.TC = M.TA;

%%% crater geometry %%%
load Erebus_res1D_geom.mat
shape = fliplr(shape);
depth = shape(1,1);

%%% load source %%%
src = load('monopole_src_1.txt');
ts = src(:,1)';
s = src(:,2)'*amp;

% formulas assume that there are an even number of time samples 
ts = ts(1:end-1); 
s = s(1:end-1);

dtintp = 0.01;
tintp = min(ts):dtintp:max(ts);
sintp = pchip(ts,s,tintp);

ts = tintp;
s = sintp;

T = max(ts); % length of time vector
N = length(ts); % number of time points
dt = T/N; % time step
Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples

%%% fourier transform source function
S = fft(s)*dt;
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
f = omega/(2*pi); % frequency vector

% convert source from mass flux to volume flux
rho = 1.22;
s_vol = s./rho;

% cross-sectional area at base of crater
A = pi*shape(1,2)^2;

% velocity at base of crater
vel = s_vol/A;

% compute source in frequency domain
S = fft(vel)*dt;

%%% resonance 1D %%%
for i = 1:length(rsta_multi)
    i
    % set parameters
    M.r = rsta_multi(i);

    % compute Green's function 
    A = resonance1d(shape, depth, freq, Nf, style, order, M);
    
    % convolve Green's function with source
    pres = A.P(1:N/2+1) .* S(1:N/2+1);
    
    % plot infrasound signal in frequency domain
    figure(5);
    subplot(num_sta,1,i); hold on;
    plot(A.f(1:N/2+1), abs(pres),'Color','k');
    %ylim([0 10]);
    xlim([0 3]);
    %ylabel('$Pa/Hz$','interpreter','latex');
    ylabel('Pa/Hz');
    
    % invert infrasound signal to time domain
    presF = [pres conj(pres(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
    presF(N/2+1) = real(presF(N/2+1)); % entry at Nyquist must be real
    presTime = ifft(presF,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain

    % plot infrasound signal in time domain
    figure(4);
    subplot(num_sta,1,i); hold on;
    plot(ts, presTime,'k');
    %xlim([0 5]);
    %ylim([-20 20]);
    %ylabel('$Pa$','interpreter','latex');
    ylabel('Pa');
    
    p_1D(i) = max(presTime);
    
    if i == num_sta
        figure(5);
        subplot(num_sta,1,i);
        xlabel('Freq (Hz)');
        
        figure(4);
        subplot(num_sta,1,i);
        xlabel('Time (s)');
    end
end