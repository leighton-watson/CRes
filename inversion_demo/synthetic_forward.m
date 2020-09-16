%% FORWARD SYNTHETIC %%
%
% Script file that performs the forward problem of calculating the
% infrasound spectra for a given crater geometry.
%
% Written by Leighton Watson
% March 12, 2020
% leightonwatson@stanford.edu // leightonmwatson@gmail.com

clear all; clc;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

path(pathdef)
addpath ../source/resonance/
addpath ../source/SBPoperators/
addpath ../source/inv/
addpath data/

save_output = 1; % logical to save output normalized amplitude spectra

%% resonance 1D properties %%

% set the parameters for the resonance1d calculations 
T = 25; % total time (s)
N = 250; % number of grid points (formulas assume even N)
dt = T/N; % time step (s)
Nyquist = 1/(2*dt); % Nyquist frequency (Hz)
Nf = N/2+1; % number of frequency samples
freq = [0 Nyquist]; % frequency range (Hz)
discrParams = [T N Nf Nyquist dt]; % save parameters into array

craterTemp = 100; % crater temperature
atmoTemp = 0; % atmospheric temperature
temp = [craterTemp, atmoTemp]; 

order = 4; % order of numerical scheme (4, 6 or 8)
style = 'baffled piston'; % acoustic radiation model ('monopole' or ' baffled piston')
M = problemParametersInv(craterTemp,atmoTemp); % problem parameters required for resonance1d

%% source properties %%
% (volumetric flow rate [m^3/s]) at base of crater
resParams = [T, N];
srcStyle = 'Gauss'; % choose from "Gauss" or "Brune"
srcAmp = 1; % amplitude of source
srcWidth = 0.2; % width of source 
[S,f,s,ts] = sourceFunction(srcAmp,srcWidth,srcStyle, resParams);

figure(1); clf;
subplot(1,2,1);
plot(ts,s,'Color','k');
xlabel('Time (s)');
ylabel('s(t)');
xlim([0 10])

subplot(1,2,2);
plot(f(1:N/2+1),abs(S(1:N/2+1))./max(abs(S(1:N/2+1))),'Color','k');
xlabel('Frequency (Hz)');
ylabel('s(\omega)');
xlim([0 3])

%% crater geometry properties %%
% load a specified geometry profile or generate new profile

%%% load specified geometry profile
load('Johnson2018');
shape = flipud(geometry);
depth = 120; % specify depth of crater

%%% generate new geometry profile
%%% first value of geomParam is depth, remaining Nx values are radial values
%%% radial values are equal spaced in depth between outlet and crater base
%%% radius is linearly interpolate on 1m intervals between specified points
% depth = 120; % specify depth of crater
% geomParam = [depth 120 80 60 50 40]; 
% shape = geomFunction(geomParam);

figure(2); clf;
plot(shape(:,2), shape(:,1),'Color','k');
set(gca,'YDir','Reverse');
hold on;
plot(-shape(:,2), shape(:,1),'Color','k');
plot([-shape(depth,2) shape(depth,2)],[depth depth],'Color','k');
ylim([0 depth]);
xlabel('Radius (m)');
ylabel('Depth (m)');

%% simulate infrasound signal %%

A = resonance1d(shape, depth, freq, Nf, style, order, M); % simulate transfer function
dP = A.P(1:N/2+1) .* S(1:N/2+1); % multiple transfer function with source

figure(3); clf;
plot(A.f(1:N/2+1), abs(dP(1:N/2+1))./max(abs(dP(1:N/2+1))),'Color','k'); % plot normalized amplitude spectra
xlim([0 3]);
xlabel('Frequency (Hz)');
ylabel('\Delta p(\omega,r)');

%% save outputs %%

if save_output == 1 
    
    dataF = A.f(1:N/2+1)';
    dataAmp = (abs(dP(1:N/2+1))./max(abs(dP(1:N/2+1))))';
    save('forwardSpectraSynthetic','dataF','dataAmp',...
        'shape','depth','srcStyle','srcAmp','srcWidth');
end