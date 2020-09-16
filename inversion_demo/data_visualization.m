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
addpath data/
addpath ../source/resonance/
addpath ../source/SBPoperators/
addpath ../source/inv/

%% load and plot inversion input %%

dataStr = 'Etna2018Phase1';
datafile = strcat(dataStr,'.mat');
load(datafile); % load data
data_freq = [dataF, dataAmp]; % format data

figure(1); clf;
subplot(1,2,2);
plot(dataF, dataAmp./max(dataAmp),'Color','k');
hold on; xlim([0 3]);
xlabel('Frequency (Hz)');
ylabel('\Delta p(\omega,r)');


%% load inversion output %%

load DataInvOut_example.mat

%% crater geometry %%    

figure(1); 

subplot(1,2,1);
%%% Initial estimate of crater geometry
shapeI = geomFunction(geomParams);
depthI = shapeI(1,1);
plot(shapeI(:,2), shapeI(:,1),'Color',cmap(1,:)); hold on;
set(gca,'YDir','Reverse'); xlabel('Radius (m)'); ylabel('Depth (m)');
plot(-shapeI(:,2), shapeI(:,1),'Color',cmap(1,:));
plot([-shapeI(1,2) shapeI(1,2)],[depthI depthI],'Color',cmap(1,:));
ylim([0 220]);
xlim([-130 130])

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

figure(1);
subplot(1,2,2);

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


%% parameter histograms %%
    
depths = [0 50 75];
param_hist = zeros(length(x_trunc),length(depths));

tic
figure(2);
for j = 1:length(x_trunc)
    disp(strcat('Percentage complete:',num2str(j/length(x_trunc)*100),'%'));
    geom_trunc = x_trunc(j,:) .* geomParams;
    shape_intp = geomFunction(geom_trunc);
    for m = 1:length(depths)
        if shape_intp(1,1) < depths(m)
            param_hist(j,m) = NaN;
        else
            param_hist(j,m) = shape_intp(end-depths(m),2);
        end
    end
end


for m = 1:length(depths)
    
    if m == 1 % plot depth
        subplot(1, length(depths)+1,1);
        h = histogram(x_trunc(:,1).*geomParams(1));
        hold on; box on;
        h.Parent.YTick = {};
        h.BinWidth = 1;
        h.FaceColor = cmap(1,:);
        h.EdgeColor = cmap(1,:);
        h.EdgeAlpha = 0;
        xlim([50 250])
    end
    
    subplot(1,length(depths)+1,m+1);
    h = histogram(param_hist(:,m));
    hold on; box on;
    h.Parent.YTick = {};
    h.BinWidth = 1;
    h.FaceColor = cmap(1,:);
    h.EdgeColor = cmap(1,:);
    h.EdgeAlpha = 0;
    
end
toc


%% parameter correlations %%
    
crater_depth = x_trunc(:,1).*geomParams(1);
mPlot = [crater_depth, param_hist];
mNames = {'Depth','Radius at z=0 m','Radius at z=50 m','Radius at z=75 m'};


PlotCorrelations(mPlot, mNames, []);


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
