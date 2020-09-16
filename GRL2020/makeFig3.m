%% FIG 3 %%
%
% Calculate peak frequency and quality factor as a function of depth for
% analytical formula (equations 2 and 4 in main text).
%
% Compute spectra for before and after onset of flank eruption. Assume a
% cylindrical crater geometry with a radius of 12 m, speed of sound of 575
% m/s and depth of 229 m (before) and 389 m (after).
%
% Written by Leighton Watson
% March 10, 2020
% leightonwatson@stanford.edu // leightonmwatson@gmail.com

clear all; clc;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

path(pathdef)
addpath data/
addpath ../source/resonance/
addpath ../source/SBPoperators/

figHand1 = figure(1); clf;
set(figHand1,'Position',[100 100 900 800]);

load EMFO_spectra.mat

%% analytical solutions for peak frequency and quality factor %%

R = [12 50 150];
L = 1:1:500;
c = 575;

% peak frequency from data
f0_before = EMFO_before.f0;
f0_after = EMFO_after.f0;

% quality factor from data
Q_before = EMFO_before.Q;
Q_after = EMFO_after.Q;

for i = 1:length(R)
    
    r = R(i);
    Leff = L + (8.*r)./(3.*pi);
    f0 = c./(4*Leff);
    Q = (4*L.^2)./(pi.*r.^2) + (32*L)./(3*pi^2.*r);
    
    figure(1); % peak frequency analytical formula
    subplot(2,2,1);
    plot(f0,L); hold on;
        
    subplot(2,2,2); % quality factor analytical formula
    semilogx(Q,L); hold on;
    
end

subplot(2,2,1);
set(gca,'YDir','Reverse');
xlim([0 1.5]);
xlabel('Peak Frequency (Hz)');
ylabel('Depth (m)');

subplot(2,2,2);
set(gca,'YDir','Reverse');
xlim([1 1e5])
xlabel('Quality Factor');

L_before = c/(4*f0_before) - (8.*R(1))./(3.*pi);
L_after = c/(4*f0_after) - (8.*R(1))./(3.*pi);

subplot(2,2,1);
vline(f0_before);
vline(f0_after);

hline(L_before);
hline(L_after);

subplot(2,2,2);
vline(Q_before);
vline(Q_after);


%% load data %%

subplot(2,2,3);
plot(EMFO_before.f, EMFO_before.spec,'k');
xlim([0 3]);
xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
hold on;

subplot(2,2,4);
plot(EMFO_after.f, EMFO_after.spec,'k');
xlim([0 3])
xlabel('Frequency (Hz)');
ylabel('Normalized Amplitude');
hold on;

%% simulate infrasound spectra with CRes %%

%%% RESONANCE 1D PARAMETERS %%%
T = 500;
N = 5000;
dt = T/N; % time step
t = [0:N-1]*dt; % vector of time grid points [s]
tshift = 0;

Nyquist = 1/(2*dt); % Nyquist frequency [Hz]
freq = [0, Nyquist]; % range of frequencies [Hz]
Nf = N/2+1; % number of frequency samples
df = Nyquist/(Nf-1) % frequency interval

%%% PROBLEM PARAMETERS %%%
M.gamma = 1.4; % ratio of heat capacities
M.r = 1000; % distance from vent to observer to calculate excess pressure at
M.R = 287; % specific gas constant for dry air

% atmospheric properties
M.TA = 550; % temperature of atmosphere [C]
M.rhoA = 1; % density [kg/m^3]
M.cA = sqrt(M.gamma*M.R*(M.TA+273.15)); % speed of sound [m/s]

% properties inside crater
M.TC = 550; % temperature inside crater [C]
M.rhoC = 1; % density [kg/m^3]
M.cC = sqrt(M.gamma*M.R*(M.TC+273.15)); % speed of sound [m/s]
M.pC = M.cC^2*M.rhoC/M.gamma; % pressure [Pa]

style = 'baffled piston'; % model of sound radiation ('baffled piston' or 'monopole')
order = 8; % order of numerical method

%%% SOURCE FUNCTION %%%
% Compute source function (volumetric flow rate [m^3/s]) at base of crater
resParams = [T, N];
srcStyle = 'Brune'; % choose from "Gauss" or "Brune"
A = 1; % source amplitude
L = 0.3; % width of source; % amplitude of source
[S,f,s,ts] = sourceFunction(A,L,srcStyle, resParams);

f0_before = EMFO_before.f0;
f0_after = EMFO_after.f0;
c = M.cC;
a = 12; % crater radius
L_before = c/(4*f0_before) - (8.*a)./(3.*pi);
L_after = c/(4*f0_after) - (8.*a)./(3.*pi);

%% CYLINDRICAL PIPE %%

%%% BEFORE ERUPTION %%%
dz = 1; % increment of depth discretization
z = 0:dz:L_before; % depth
depth = max(z);
r = a*ones(size(z)); % vector of radius
geometry = [z' r']; % matrix of geometry
shape = flipud(geometry); % reformat geometry matrix to be read by solver


A = resonance1d(shape, depth, freq, Nf, style, order, M); % solve transfer function
transFunc = A.P;
dp = transFunc .* S(1:N/2+1);

figure(1); subplot(2,2,3);

plot(A.f, abs(dp)./max(abs(dp)));
xlim([0 3]);

[f0_simb,Q_simb] = resPeakProps(A.f, abs(dp)./max(abs(dp)));

%%% AFTER ERUPTION %%%
dz = 1; % increment of depth discretization
z = 0:dz:L_after; % depth
depth = max(z);
r = a*ones(size(z)); % vector of radius
geometry = [z' r']; % matrix of geometry
shape = flipud(geometry); % reformat geometry matrix to be read by solver

A = resonance1d(shape, depth, freq, Nf, style, order, M); % solve transfer function
transFunc = A.P;
dp = transFunc .* S(1:N/2+1);

figure(1); subplot(2,2,4);

plot(A.f, abs(dp)./max(abs(dp)));
xlim([0 3]);

[f0_sima,Q_sima] = resPeakProps(A.f, abs(dp)./max(abs(dp)));

