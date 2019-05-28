%% MAKE FIG %%
%
% Make figure showing geometry and infrasound time series for several
% different slope angles away from crater rim

clear all; clc; %close all;

%add paths
addpath ../source/SBPoperators/
addpath ../source/resonance/
addpath Data/Fig9_RadiationConicalVolcanoes/


set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',20);
cmap = get(gca,'ColorOrder');

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 1300 700]);
%% INFRA FDTD %%
% 
% Plot infraFDTD simulations for several slope angles away from the crater
% rim. Plot geometry and infrasound time series.

%path_str = {'geom000','geom015','geom030','geom060'}; % directory
path_str = {'geom000','geom015','geom030'};
path_num = length(path_str);

nx = 201; % number of grid points in x direction
ny = 201; % number of grid points in y direction
%nz = [83 83 84 128]; % number of grid points in z direction

% infrasound stations
sta_str1 = 'ST';
sta_str2 = '_pressure.txt';
num_sta = 4;
idx_plot = [2 3 5 6]; % index for plotting 

% receiver position
%y_pos = [1000 851.4 739.2];
y_pos = [1005 925 855];
r_horz = 900;

% determine volcano apex
thetad = [0 15 30];
crater_apex = tand(thetad)*100;

for i = 1:length(path_str)
    
    %%% add path %%%
    path(pathdef); % clear previous paths
    path_str_load = strcat('Data/Fig9_RadiationConicalVolcanoes/',path_str{i},'/'); % form path string
    addpath(path_str_load); % add path
    
    %%% geometry %%%
    elev = load('sur_coords.txt');
    x = elev(:,1);
    y = elev(:,2);
    z = elev(:,3);
    
    X = reshape(x,nx,ny);
    Y = reshape(y,nx,ny);
    Z = reshape(z,nx,ny);
    
    % plot profile 
    idx = ceil(nx/2);
    figure(1); 
    subplot(path_num,3,3*i-2);
    h = plot(Y(idx,:), Z(idx,:),'k');
    ylim([750 1150])
    set(h.Parent,'YTick',[800 900 1000 1100]);
    set(h.Parent,'YTickLabel',{'200','','0',''});
    hold on;
    xlabel('Radius (m)');
    ylabel('Depth (m)');
    %grid on
    plot(r_horz, y_pos(i)+20,'v','MarkerFaceColor','k',...
        'MarkerEdgeColor','k','MarkerSize',14,'LineWidth',1);
    
    % plot apex of crater
    plot(500,1000+crater_apex(i),'o','MarkerFaceColor','k',...
        'MarkerEdgeColor','k','MarkerSize',12,'LineWidth',1);
    
    %%% time series %%%
    for j = 4; %1:num_sta
    
        % load data
        sta_str = strcat(sta_str1,num2str(j),sta_str2);
        sta = load(sta_str);
        
        % plot data
        figure(1); 
        subplot(path_num,3,[3*i-1 3*i]);
        %plot(sta(:,1), sta(:,2),'Color',cmap(i,:));
        h = plot(sta(:,1), sta(:,2),'k');
        set(h.Parent,'YTick',[-2 -1 0 1 2]);
        set(h.Parent,'YTickLabel',{'-2','','0','','2'});
        set(h.Parent,'XTick',[0:1:10]);
        set(h.Parent,'XTickLabel',{'0','','2','','4','','6','','8','','10'});
        xlabel('Time (s)');
        ylabel('$\Delta p(t,r)$','interpreter','latex');
        xlim([0 10]);
        hold on;
        %grid on
        ylim([-2 2])
        
        
    end
        
    
end




%% RES 1D %%
%
% Compare infraFDTD simulations to res1D simulations for radiation into 
% a half space. Can propose analytical corrections for the reduced
% amplitude when the volcano slopes downwards from the crater.

addpath ../source/SBPoperators/
addpath ../source/resonance/
addpath Data/Fig9_RadiationConicalVolcanoes/
addpath Functions/

%r = [100 200 300 400]; % radial distance to receiver
r_horiz = 400;
%thetad = [0 15 30 60];
thetad = [0 15 30];
theta = thetad*(2*pi/360);

%y = 1000 - [1000 920 851.4 739.2]; %r_horiz*tan(theta);

%y = 1000 - [1000 851.4 739.2];
%r = (r_horiz.^2 + y.^2).^(0.5);
r = 400./cosd(thetad);


gamma = 1; % isothermal to match infraFDTD solutions
rhoA = 1.22; % atmospheric density
rhoC = 1.22; % crater air density
cA = 340; % atmospheric speed of sound
cC = 340; % crater air speed of sound
%style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
style = 'monopole';
order = 8; % order of numerical method

% load crater geometry
shape = load('res1d_geometry.txt');
depth = shape(1,1);

% load source function
src = load('monopole_src_1.txt');
ts = src(:,1)';
s = src(:,2)';

% formulas assume that there are an even number of time samples
ts = ts(1:end-1); 
s = s(1:end-1);

T = max(ts); % length of time vector
N = length(ts); % number of time points
dt = T/N; % time step
Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples
omega = [0:N/2 (-N/2+1):-1]/N*2*pi/dt; % angular frequency vector
f = omega/(2*pi); % frequency vector

% convert from mass flux to volume flux
s_vol = s/rhoC;

% cross-sectional area at base of crater
A = pi*shape(1,2)^2;

% convert to velocity at base of crater
s_vel = s_vol/A;

% Fourier transform source function
S = fft(s_vel)*dt;

%%% Compute Green's function and synthetic pressure. Iterate over the
%%% different distances

for k = 1:path_num %length(r)
    
    % load problem parameters
    M = problemParameters_validation(gamma, r(k)); % load problem parameters
    M.rhoA = rhoA;
    M.rhoC = rhoC;
    M.cA = cA;
    M.cC = cC;
    M.pC = M.cC^2*M.rhoC/M.gamma;
    M.TA = M.cA.^2/(M.gamma*M.R) - 273.15;
    M.TC = M.TA;
    
    % compute Green's function
    A = resonance1d_slope(shape, depth, freq, Nf, style, order, M, thetad(k));
    
    % convolve Green's function with source
    pres = A.P(1:N/2+1) .* S(1:N/2+1);
    
    % invert infrasound signal to the time domain
    presF = [pres conj(pres(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
    presF(N/2+1) = real(presF(N/2+1)); % entry at Nyquist must be real
    presTime = ifft(presF,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain
    
    % plot infrasound signal in time domain
    figure(1);
    subplot(path_num,3,[3*k-1 3*k]);
    %plot(ts, presTime,'k');
    %plot(ts+r(k)/M.cA, presTime,'Color',cmap(k,:),'LineStyle',':');
    h = plot(ts+r(k)/M.cA, presTime,'Color',cmap(k,:),'LineStyle','-');
    h.Color(4) = 1;
    
    h = plot([0 1.5],[0 0],'Color',cmap(k,:),'LineStyle','-');
    h.Color(4) = 1;
    
    %%% compare to no correction to analytical radiation model %%%
    M.r = r_horiz;
    A = resonance1d_slope(shape, depth, freq, Nf, style, order, M, 0);
    pres = A.P(1:N/2+1) .* S(1:N/2+1);
    presF = [pres conj(pres(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
    presF(N/2+1) = real(presF(N/2+1)); % entry at Nyquist must be real
    presTime = ifft(presF,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain
    subplot(path_num,3,[3*k-1 3*k]);
    h = plot(ts+r_horiz/M.cA, presTime,'Color',cmap(k,:),'LineStyle',':');
    h.Color(4) = 0.6;

   
    
end
