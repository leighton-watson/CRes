%% MAKE FIG %%
% IMPEDANCE CONTRAST
%
% Make plot of impedance contrast and example infrasound signals to include
% in the appendix of revised version of Infrasonic resonance of volcanic
% craters paper

% add paths
addpath ../source/SBPoperators/
addpath ../source/resonance/

clear all; clc; %close all;
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',20);
cmap = get(gca,'ColorOrder');

figHand1 = figure(1); %clf;
figHand2 = figure(2); clf;
figHand3 = figure(3); clf;

set(figHand1,'Position',[100 100 600 700]);
set(figHand2,'Position',[1000 100 600 300]);
set(figHand3,'Position',[1000 400 600 300]);

% parameters
gamma = 1.4; % ratio of heat capacities
R = 287; % specific gas constant for dry air

% crater geometry at outlet
a = 100; % radius of crater outlet [m]
A = pi*a^2; % area of crater outlet [m^2]

% atmospheric properties
Ta = 20; % temperature of atmosphere [c]
rhoa = 1; % density of atmosphere [kg/m^3]
ca = sqrt(gamma*R*(Ta+273.15)); % speed of sound [m/s]
pa = ca^2*rhoa/gamma; % pressure [Pa]

% frequency
freq = [0 10]; % frequency range [Hz]
Nf = 1000; % number of frequency samples
f = linspace(freq(1),freq(end),Nf); % frequency vector [Hz]
k = 2*pi*f/ca; % wave number using speed of sound in atmosphere

% properties at crater outlet
Tc = [0]; % temperature inside crater [c]
Tk = Tc+273.15; % temperature inside crater in Kelvin
rhoc = pa/(R*Tk); % density inside crater [kg/m^3]
cc = sqrt(gamma*R*Tk); % speed of sound [m/s]
ZN = rhoc*cc/A;

for j = 1:length(f)
    ZT(j) = flanged_opening(f(j),rhoa,ca,a);
end

figure(1); clf;
semilogy(k.*a, abs(ZT)/abs(ZN));
hold on;

xlabel('$ka$','interpreter','latex')
ylabel('$|Z_T|/|Z_N|$','interpreter','latex')
hline(1)
xlim([0.02 2])

%% simulate infrasound signal
 
T = 100; % total time [s]
N = 2000; % number of grid points in time (formulas below assume even N)
dt = T/N; % time step
t = [0:N-1]*dt; % vector of time grid points [s]

Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples

M = problemParameters(); % property values are defined in resonance/problemParameters.m

M.gamma = gamma; % ratio of heat capacities
M.r = 1000; % distance from vent to observer to calculate excess pressure at
M.R = R; % specific gas constant for dry air

% atmospheric properties
M.TA = Ta; % temperature of atmosphere
M.rhoA = rhoa; % density [kg/m^3]
M.cA = sqrt(M.gamma*M.R*(M.TA+273.15)); % speed of sound [m/s]
 
% properties inside crater
M.TC = Tc; % temperature inside crater
M.rhoC = rhoc; % density [kg/m^3]
M.cC = sqrt(M.gamma*M.R*(M.TC+273.15)); % speed of sound [m/s]
M.pC = M.cC^2*M.rhoC/M.gamma; % pressure [Pa]

style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method

%%% Source Model %%%
amp = 1;
L = 0.4;
srcStyle = 'Gauss';
resParams = [T N];
[S, f, s, ts] = sourceFunction(amp, L, srcStyle, resParams);

%%% Crater Properties %%%
radius_factor = [1 2];
depth = 200;

for i = 1:length(radius_factor)
   
    shape1 = [1000:-1:0]';
    shape2 = 200*ones(size(shape1));
    shape = [shape1 shape2];
    shape(:,2) = shape(:,2).*radius_factor(i);
    
    rhoc = M.rhoC; % density inside crater [kg/m^3]
    cc = M.cC; % speed of sound [m/s]
    area = pi*shape(end,2)^2;
    %ZN = rhoc*cc/area;
    ZN = M.rhoA*M.cA/area;
    
    for j = 1:length(f)
       ZT(j) = flanged_opening(f(j),rhoa,ca,a);        
    end
    
        Zratio(i) = abs(ZT(end))/abs(ZN)
        figure(1); hline(Zratio(i));
        A = resonance1d(shape, depth, freq, Nf, style, order, M);
        
        pres = A.P(1:N/2+1) .* S(1:N/2+1);
        [f0(i), Q(i)] = resPeakProps(A.f, abs(pres)./max(abs(pres)))
        
        %ka(i) = 2*pi*f0(i)/M.cC*shape2(1);
        %vline(ka(i))
        
       
        
        figure(1+i); 
        subplot(1,2,1);
        hold on;
        h = plot(A.f, abs(pres)./max(abs(pres)),'Color',cmap(i+1,:)); % plot vs frequency
        set(h.Parent,'YTick',[]);
        
        %plot(ka, abs(pres)./max(abs(pres)));
        xlim([0 1])
        box on
        xlabel('Frequency (Hz)');
        ylabel('$\Delta p(\omega,r)$','interpreter','latex')
        
        % plot in time domain
        sig_pos = pres;
        sig_full = [sig_pos conj(sig_pos(end-1:-1:2))]; % reflect about f = 0. Take complex conjugate for negative frequencies
        sig_full(N/2+1) = real(sig_full(N/2+1)); % entry at Nyquist must be real
        sig_time = ifft(sig_full,'symmetric')/dt; % inverse Fourier transform to obtain Greens function in time domain

    
        figure(1+i); 
        subplot(1,2,2); hold on;
        h = plot(t-20, sig_time./max(abs(sig_time)),'Color',cmap(i+1,:));
        xlim([0 40])
        set(h.Parent,'YTick',[]);
        box on
        xlabel('Time (s)');
        ylabel('$\Delta p(t,r)$','interpreter','latex')
end

