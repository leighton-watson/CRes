%% INVERSION ETNA %%
%
% Script file that performs MCMC inversion to invert harmonic infrasound
% observations for crater geometry using real data from Mount Etna
%
% Written by Leighton Watson
% March 10, 2020
% leightonwatson@stanford.edu // leightonmwatson@gmail.com

clear all; clc;

path(pathdef)
addpath ../../source/resonance/
addpath ../../source/SBPoperators/
addpath ../../source/inv/
addpath ../data/


%% resonance 1D %%

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

%% data %%


% filter properties
filterband = [0.15 4.8]; % frequency band to filter
filterorder = 2; % order of butterworth filter
Fs = 10; % sampling frequency
filterProps = [filterband, filterorder, Fs]; % filter properties - same as for data

% data
dataStr = 'EMFO_before';
datafile = strcat(dataStr,'.mat');
load(datafile); % load data
data_freq = [dataF', dataAmp']; % format data

%% inversion parameters %%

nIter = 100; % number of steps

dx = 0.05; % step size % use step size of 0.05 for paper inversions

freqLim = 2; % high cut frequency limit for misfit function (Hz)

%%% geometry parameters %%%
geomFlag = 1; % invert for geometry (boolean, 0 = no, 1 = yes)
geomR0 = 100; % radius of initial cylinder
geomDepth = 200; % depth 
geomParams = [geomDepth geomR0 geomR0 geomR0 geomR0 geomR0]; % first value is depth, other values are radius points that are equally spaced
geomLowerBnds = [50 80 1 1 1 1];
geomUpperBnds = [300 120 120 120 120 120];
nx = length(geomParams)-1; % number of geometry parameters

%%% source parameters %%%
srcFlag = 0; %invert for source (boolean, 0 = no, 1 = yes)
srcParams = 0.3;
srcStyle = 'Brune';
srcUpperBnds = [5];
srcLowerBnds = [0.01];

%%% format parameters %%%
if geomFlag && srcFlag
    upperBnds = [geomUpperBnds srcUpperBnds];
    lowerBnds = [geomLowerBnds srcLowerBnds];
    params = [geomParams srcParams];
    
elseif geomFlag
    upperBnds = geomUpperBnds;
    lowerBnds = geomLowerBnds;
    params = geomParams;
    
elseif srcFlag
    upperBnds = srcUpperBnds;
    lowerBnds = srcLowerBnds;
    params = srcParams;
   
end

%% mcmc inversion %%

tic % start timing
[x, misfit, simSpec, f, count] = mcmc_spec(nIter, dx, ... % perform MCMC inversion
    geomParams, geomFlag, srcParams, srcFlag, srcStyle,...
    lowerBnds, upperBnds, discrParams, temp, ...
    filterProps, data_freq, freqLim);
toc % finish timing
disp(['Elapsed time is ',num2str(toc/60), ' minutes.']); % display timing in minutes


%% save outputs %%

% format spectra
spec.int = simSpec(1,:); % initial spectra
spec.fin = simSpec(end,:); % final spectra
burn_in = ceil(count/10); % remove the first 10% of successful samples to remove effect of initial conditions
simSpec_trunc = simSpec(burn_in:end,:); % remove burn in to reduce sensitivity to initial conditions
simSpec_mean = mean(abs(simSpec_trunc)); % mean spectra
spec.mean = abs(simSpec_mean)./max(abs(simSpec_mean)); % normalize spectra

pathname = strcat(pwd,'/invOutput'); % directory to save outputs
filename = strcat('InvOut','_',dataStr,'_',... % file name
    'Nit',num2str(nIter),'_',...
    'R0',num2str(geomR0),'m_','D',num2str(geomDepth),'m_',...
    'T',num2str(craterTemp),'C_',...
    'freqLim',num2str(freqLim),'_nx',num2str(nx));
matfile = fullfile(pathname,filename); % path and file for saving outputs

save(matfile,'x','misfit','spec','f','count',...
    'params','geomParams','srcParams','geomFlag','srcFlag','srcStyle',...
    'lowerBnds','upperBnds','nIter','dx','freqLim',...
    'discrParams','M','filterProps');



