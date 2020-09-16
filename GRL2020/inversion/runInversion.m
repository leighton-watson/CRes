%% INVERSION ETNA %%
%
% Script file that performs MCMC inversion to invert harmonic infrasound
% observations for crater geometry using real data from Mount Etna
%
% Written by Leighton Watson
% March 10, 2020
% leightonwatson@stanford.edu // leightonmwatson@gmail.com

clear all; clc;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

path(pathdef)
addpath source/resonance/
addpath source/SBPoperators/
addpath source/inv/

for i = 1:2
    figure(i); clf;
end

save_output = 1; % logical that determines if outputs are saved or not
plot_output = 1; % logical that determines if outputs are plotted

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

% load EMFO_data.mat;
% dataStr = 'before';
% 
% if strcmp(dataStr,'before')
%    
%     data.t_hr = EMFO_before.t; % time in hours
%     data.t = data.t_hr*3600; % time in seconds
%     data.p = EMFO_before.p; % pressure
%     
% elseif strcmp(dataStr,'after')
%     
%     data.t_hr = EMFO_after.t; % time in hours
%     data.t = data.t_hr*3600; % time in seconds
%     data.p = EMFO_after.p; % pressure
%     
% else
%     
%     error('Incorrect string: select from ''before'' or ''after''.');
%     
% end
% 
% 
% filter properties
filterband = [0.15 4.8]; % frequency band to filter
filterorder = 2; % order of butterworth filter
Fs = 10; % sampling frequency
filterProps = [filterband, filterorder, Fs]; % filter properties - same as for data

% figure(1); clf;
% plot(data.t, data.p);
% hold on;
% 
% L = length(data.p);
% NFFT = 2^nextpow2(L);
% data.F = fft(data.p,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% 
% figure(2); clf;
% plot(f, abs(data.F(1:NFFT/2+1)));
% 
% data.pfilt = bandpass_butterworth(data.t,filterband,Fs,filterorder);
% 
% figure(1); hold on;
% plot(data.t, data.pfilt);
% 
% 
% figure(2); hold on;
% L = length(data.pfilt);
% NFFT = 2^nextpow2(L);
% data.F = fft(data.pfilt,NFFT)/L;
% f = Fs/2*linspace(0,1,NFFT/2+1);
% plot(f, abs(data.F(1:NFFT/2+1)));
% 
% return



% 
dataStr = 'EMFO_before';
datafile = strcat(dataStr,'.mat');
load(datafile); % load data
data_freq = [dataF', dataAmp']; % format data

%%%% NEED TO CHANGE FILTER PARAMETERS %%%


%% inversion parameters %%

nIter = 100; % number of steps

dx = 0.1; % step size % use step size of 0.05 for paper inversions

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

if save_output == 1
    save(matfile,'x','misfit','spec','f','count',...
    'params','geomParams','srcParams','geomFlag','srcFlag','srcStyle',...
    'lowerBnds','upperBnds','nIter','dx','freqLim',...
    'discrParams','M','filterProps');
end

%% plot outputs %%

if plot_output == 1
    
    %%% crater geometry %%%
    
    figure(1);
    
    %%% Initial estimate of crater geometry
    shapeI = geomFunction(geomParams);
    depthI = shapeI(1,1);
    plot(shapeI(:,2), shapeI(:,1),'Color',cmap(1,:)); hold on;
    set(gca,'YDir','Reverse'); xlabel('Radius (m)'); ylabel('Depth (m)');
    plot(-shapeI(:,2), shapeI(:,1),'Color',cmap(1,:));
    plot([-shapeI(1,2) shapeI(1,2)],[depthI depthI],'Color',cmap(1,:));
    
    %%% Mean estimate of crater geometry
    burn_in = ceil(count/10); % remove the first 10% of successful samples to remove effect of initial conditions
    x_out = x(:,1:length(geomParams)); % extract output geometry parameters
    x_trunc = x_out(burn_in:end,:); % remove burn-in values
    x_mean = mean(x_trunc); % average parameters
    geomParams_mean = x_mean.*geomParams; % convert to physical values
    shapeF = geomFunction(geomParams_mean);
    depthF = shapeF(1,1);
    plot(shapeF(:,2), shapeF(:,1),'Color',cmap(2,:));
    plot(-shapeF(:,2), shapeF(:,1),'Color',cmap(2,:));
    plot([-shapeF(1,2) shapeF(1,2)],[depthF depthF],'Color',cmap(2,:));

    %%% spectra %%%

    figure(2);
    
    %%% Data
    plot(dataF, abs(dataAmp)./max(abs(dataAmp)),'k');
    hold on; xlim([0 3]);
    xlabel('Frequency (Hz)');
    ylabel('Normalized Amplitude Spectra');
    
    %%% Initial estimate of spectra
    plot(f, abs(spec.int),'Color',cmap(1,:));
    
    %%% Final estimate of spectra
    plot(f, abs(spec.fin),'Color',cmap(2,:));
    
    %%% Mean estimate of spectra
    plot(f, abs(spec.mean),'Color',cmap(3,:));
    
    %%% Spectra for mean estimate of crater geometry
    % source
    [S, ~, ~, ~] = sourceFunction(1, srcParams, srcStyle, discrParams); % compute source in frequency domain
    % resonance 1d properties
    style = 'baffled piston'; % acoustic radiation model ('monopole' or ' baffled piston')
    order = 4; % order of numerical scheme (4, 6 or 8)
    N = discrParams(2); % number of grid points (in time)
    Nf = discrParams(3); % number of frequency samples
    Nyquist = discrParams(4); % Nyquist frequency
    dt = discrParams(5); % time step
    freq = [0 Nyquist]; % frequency vector
    % resonance 1d
    res = resonance1d(shapeF, depthF, freq, Nf, style, order, M); % compute transfer function
    sim.f = res.f; % frequency vector
    sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1)); % convolve transfer function and source function and normalize amplitudes
    % convert to time domain
    sim.pFreq = [sim.P conj(sim.P(end-1:-1:2))];
    sim.pFreq(N/2+1) = real(sim.pFreq(N/2+1));
    sim.pTime = ifft(sim.pFreq,'symmetric')/dt;
    % apply filtering
    filterband = [filterProps(1) filterProps(2)]; % frequency band to filter
    filterorder = filterProps(3); % order of butterworth filter
    Fs = filterProps(4);
    sim.pTimeFilt = bandpass_butterworth(sim.pTime,filterband,Fs,filterorder);
    % convert back to frequency domain
    L = length(sim.pTimeFilt);
    NFFT = 2^nextpow2(L);
    sim.pFreqFilt = fft(sim.pTimeFilt,NFFT)/L;
    sim.pFreqAbsNorm = abs(sim.pFreqFilt(1:NFFT/2+1))./max(abs(sim.pFreqFilt(1:NFFT/2+1)));
    sim.fFilt = Fs/2*linspace(0,1,NFFT/2+1);
    % plot
    plot(sim.fFilt(1:N/2+1), abs(sim.pFreqAbsNorm(1:N/2+1)),'Color',cmap(4,:));
    legend('Data','Initial','Final','Mean Spectra','Mean Geometry');

%     %%% misfit %%%
%     figure(3);
%     plot(misfit); hold on;
%     ylabel('Misfit');
%     xlabel('Number of Successful Iterations');
%     vline(burn_in);
%     xlim([1 count])
% 
%     %%% parameter histograms %%%
%     figure(4);
%     x_trunc = x(burn_in:end,:); % remove burn-in values
%     k = length(params);
%     kGeom = length(geomParams);
%     kSrc = length(srcParams);
%     for i = 1:k
%         
%         paramPlot = params(i)*x_trunc(:,i);
%         subplot(1,k,i); hold on; box on;
%         histogram(paramPlot);
%         vline(lowerBnds(i));
%         vline(upperBnds(i));
%         xlim([lowerBnds(i) upperBnds(i)])
%     end
% 
% 	%%% parameter time evolution %%%
%     figure(5);
%     for i = 1:k
%         paramPlot = params(i)*x_trunc(:,i);
%         subplot(1,k,i); hold on; box on;
%         plot(paramPlot);
%         vline(burn_in);
%         xlim([1 count])
%     end
    
    
end
