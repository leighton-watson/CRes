%% INVERSION VISUALIZATION %%
%
% Script file that visualizes outputs from inversion scheme
%
% Written by Leighton Watson
% March 12, 2020
% leightonwatson@stanford.edu // leightonmwatson@gmail.com

clear all; clc;
cmap = get(gca,'ColorOrder');
set(0,'DefaultLineLineWidth',3);
set(0,'DefaultAxesFontSize',18);

path(pathdef)
addpath invOutput/
addpath ../source/resonance/
addpath ../source/SBPoperators/
addpath ../source/inv/

%% load and plot inversion input %%

dataStr = 'forwardSpectraSynthetic';
datafile = strcat(dataStr,'.mat');
load(datafile); % load data
data_freq = [dataF, dataAmp]; % format data

figure(1); clf;
didx = find(shape(:,1) == depth);
plot(shape(end:-1:didx,2), shape(end:-1:didx,1),'Color','k');
set(gca,'YDir','Reverse');
hold on;
plot(-shape(end:-1:didx,2), shape(end:-1:didx,1),'Color','k');
plot([-shape(depth,2) shape(depth,2)],[depth depth],'Color','k');
xlabel('Radius (m)');
ylabel('Depth (m)');

figure(2); clf;
plot(dataF, dataAmp,'Color','k');
hold on; xlim([0 3]);
xlabel('Frequency (Hz)');
ylabel('\Delta p(\omega,r)');


%% load inversion output %%

load SyntheticInvOut_Nit1000_R0100m_D150m_T100C_freqLim3_nx7;

%% crater geometry %%    

figure(1); 

%%% Initial estimate of crater geometry
shapeI = geomFunction(geomParams);
depthI = shapeI(1,1);
plot(shapeI(:,2), shapeI(:,1),'Color',cmap(1,:)); hold on;
set(gca,'YDir','Reverse'); xlabel('Radius (m)'); ylabel('Depth (m)');
plot(-shapeI(:,2), shapeI(:,1),'Color',cmap(1,:));
plot([-shapeI(1,2) shapeI(1,2)],[depthI depthI],'Color',cmap(1,:));
ylim([0 170]);
xlim([-120 120])

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

%%% Standard deviation 
%standard deviation of geometry (positive)
xStd = std(x_trunc); % compute standard deviation once burn in has been removed
geomStd1 = (x_mean + xStd) .* geomParams; % convert to physical values
shapeStd1 = geomFunction(geomStd1);
depthStd1 = shapeStd1(1,1);

plot(shapeStd1(:,2), shapeStd1(:,1),...
    'Color', cmap(2,:),'LineStyle',':'); hold on;
plot(-shapeStd1(:,2), shapeStd1(:,1),...
    'Color', cmap(2,:),'LineStyle',':');
plot([-shapeStd1(1,2) shapeStd1(1,2)],[depthStd1 depthStd1],...
    'Color',cmap(2,:),'LineStyle',':');

%standard deviation of geometry (negative)
geomStd2 = (x_mean - xStd) .* geomParams; % convert to physical values
shapeStd2 = geomFunction(geomStd2);
depthStd2 = shapeStd2(1,1);

plot(shapeStd2(:,2), shapeStd2(:,1),...
    'Color', cmap(2,:),'LineStyle',':'); hold on;
plot(-shapeStd2(:,2), shapeStd2(:,1),...
    'Color', cmap(2,:),'LineStyle',':');
plot([-shapeStd2(1,2) shapeStd2(1,2)],[depthStd2 depthStd2],...
    'Color',cmap(2,:),'LineStyle',':');


        
%% spectra %%

figure(2);

%%% Initial estimate of spectra
plot(f, abs(spec.int),'Color',cmap(1,:));
hold on;

%%% Final estimate of spectra
plot(f, abs(spec.fin),'Color',cmap(2,:));

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
% plot
plot(sim.f(1:N/2+1), abs(sim.P(1:N/2+1)),'Color',cmap(3,:));
legend('Data','Initial','Final','Mean Spectra','Mean Geometry');

%% misfit %%
figure(3);
plot(misfit); hold on;
ylabel('Misfit');
xlabel('Number of Successful Iterations');
vline(burn_in);
xlim([1 count])

%% parameter histograms %%
figure(4);
x_trunc = x(burn_in:end,:); % remove burn-in values
k = length(params);
kGeom = length(geomParams);
kSrc = length(srcParams);
for i = 1:k
    
    paramPlot = params(i)*x_trunc(:,i);
    subplot(1,k,i); hold on; box on;
    histogram(paramPlot);
    vline(lowerBnds(i));
    vline(upperBnds(i));
    xlim([lowerBnds(i) upperBnds(i)])
end

%% parameter time evolution %%
figure(5);
for i = 1:k
    paramPlot = params(i)*x_trunc(:,i);
    subplot(1,k,i); hold on; box on;
    plot(paramPlot);
    vline(burn_in);
    xlim([1 count])
end
