%% MCMC EVAL %%
%
% Use MCMC algorithm to invert for volcano properties (crater geometry and
% source mechanism)

clear all; 
clc;
path(pathdef);
set(0,'DefaultLineLineWidth',2);
set(0,'DefaultAxesFontSize',16);
cmap = get(gca,'ColorOrder');

path(pathdef)
addpath ../../source/resonance/
addpath ../../source/SBPoperators/
%addpath resonance_test
%addpath SBPoperators


%% DATA %%

load('data_late.mat');
figure(1); clf;
plot(data_freq(:,1), abs(data_freq(:,2))./max(abs(data_freq(:,2))),'k');
xlabel('Frequency (Hz)');
xlim([0 5])

figure(2); clf;
plot(data_time(:,1), data_time(:,2)./max(data_time(:,2)),'k');
hold on;


%% Resonance Initialize %%
%
% Set the parameters for the resonance1d calculations (frequency range,
% number of samples etc.)

T = 25; % total time (s)
N = 250; % number of grid points (formulas assume even N)
dt = T/N; % time step (s)
Nyquist = 1/(2*dt); % Nyquist frequency (Hz)
Nf = N/2+1; % number of frequency samples

% save parameters into array
discrParams = [T N Nf Nyquist];

%% Parameter Initialization %%
%
% Initialize the variables to be calculated with the inverse scheme

dx = 0.01;

%%% misfit parameters %%%
%
% the misfit function is:
% misfit = delta*abs(f0_data - f0_sim)/f0_data + gamma*abs(Q_data - beta*Q_sim)/Q_data
beta = 1; % accounts for less attenuation in simulation than data and how 
delta = 3; % weighting of resonant frequency term
gamma = 1; % weighting of quality factor
misfitParams = [beta delta gamma];

%%% geometry parameters %%%
geomFlag = 1; % invert for geometry (boolean, 0 = no, 1 = yes)
%geomParams = [150 107 42 30 24 16];
%geomParams = [80 100 70 50 30 30];
geomParams = [150 100 100 100 100 100];
geomLowerBnds = geomParams*0.1;
geomUpperBnds = geomParams*2;
%geomLowerBnds = [50 100 30 20 10 10];
%geomUpperBnds = [200 115 100 80 50 50];
%geomUpperBnds = [200 180 120 90 80 80];
geomStyle = 'intp';

%%% source parameters %%%
srcFlag = 0; %invert for source (boolean, 0 = no, 1 = yes)
srcParams = 0.3;
srcStyle = 'Brune';
srcUpperBnds = [0.5];
srcLowerBnds = [0.1];

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


%% MCMC INVERSION %%

nIter = 1000; % number of steps


tic % start timing
[x, f0, Q, misfit, simSpec, f, count] = mcmc(nIter, dx, ... % perform MCMC inversion
    geomParams,geomFlag, geomStyle, srcParams, srcFlag, srcStyle,...
    lowerBnds,upperBnds, discrParams, misfitParams, ...
    data_freq);
toc % finish timing

%% VISUALIZATION %%

%%% resonant frequency and quality factor %%%
figure(3); clf;

% compute resonant frequency and quality factor
dataSpec.P = data_freq(:,2);
dataSpec.f = data_freq(:,1);
[f0d Qd] = resPeakPropsInversion(dataSpec,'data');

% resonant frequency vs iteration number
subplot(3,3,[1 2]);
plot(f0);
ylabel('Resonant Freq (Hz)');
xlabel('Iteration Number');
hold on;
hline(f0d);

% resonant frequency histogram
subplot(3,3,3);
histogram(f0);
xlabel('Resonant Freq (Hz)');
hold on;
vline(f0d);

% quality factor vs iteration number
subplot(3,3,[4 5]);
plot(Q);
ylabel('Quality Factor');
xlabel('Iteration Number');
hold on;
hline(Qd/beta);

% quality factor histogram
subplot(3,3,6);
histogram(Q);
hold on;
vline(Qd/beta);

% misfit vs iteration number
subplot(3,3,[7 8]);
plot(misfit);
xlabel('Iteration Number');
ylabel('Misfit');

%%% inversion outputs %%%
figure(4); clf;
geomLength = length(geomParams);
xGeom = x(:,1:geomLength);
paramOutput = xGeom .* geomParams;
for i = 1:geomLength
   subplot(1,geomLength,i);
   histogram(paramOutput(:,i))
   xlim([geomLowerBnds(i) geomUpperBnds(i)]);
end

%%% source %%%
if srcFlag
    figure(7); clf;
    plot(x(:,end)*srcParams);
end
hline(srcParams);


%%% crater geometry %%%
figure(5); clf;
[m,~] = size(x);
alpha_plot = 0.7;
for i = m:-1:1
    shape = geomFunction(x(i,1:geomLength).*geomParams,geomStyle);
    depth = shape(1,1);
    
    plot(shape(:,2), shape(:,1),'Color',[1 1 1]*alpha_plot);
    set(gca,'YDir','Reverse'); hold on;
    plot(-shape(:,2), shape(:,1),'Color',[1 1 1]*alpha_plot);
    plot([-shape(1,2) shape(1,2)],[depth depth],'Color',[1 1 1]*alpha_plot)

end

% plot initial crater geometry
shapeI = geomFunction(geomParams,geomStyle);
depthI = shape(1,1);
plot(shapeI(:,2), shapeI(:,1),'Color',cmap(2,:));
plot(-shapeI(:,2), shapeI(:,1),'Color',cmap(2,:));
plot([-shapeI(1,2) shapeI(1,2)],[depthI depthI],'Color',cmap(2,:));

% plot average crater geometry after truncating burn-in values
burn_in = ceil(count/10); % remove the first 10% of successful samples to remove effect of initial conditions
x_out = x(:,1:length(geomParams)); % extract output geometry parameters
x_trunc = x_out(burn_in:end,:); % remove burn-in values
x_mean = mean(x_trunc); % average parameters
geomParams_mean = x_mean.*geomParams; % convert to physical values

% compute crater geometry for mean output valiues
shapeF = geomFunction(geomParams_mean,geomStyle);
depthF = shapeF(1,1);
figure(5);
plot(shapeF(:,2), shapeF(:,1),'Color',cmap(1,:));
plot(-shapeF(:,2), shapeF(:,1),'Color',cmap(1,:));
plot([-shapeF(1,2) shapeF(1,2)],[depthF depthF],'Color',cmap(1,:));

%%% SfM crater geometry from Villarrica paper
SfM = load('SfM_Geom.csv');
for i = 1:4
    plot(SfM(:,i+1),SfM(:,1),'Color','k','LineStyle','--');
    plot(-SfM(:,i+1),SfM(:,1),'Color','k','LineStyle','--');
end

%%% spectra %%%
figure(6); clf;
for i = 1:m
    plot(f, abs(simSpec(i,:)),'Color',[1 1 1]*alpha_plot);
    hold on;
end

% data
plot(data_freq(:,1), data_freq(:,2),'Color','k');
xlim([0 5]);

% initial estimate
plot(f, abs(simSpec(1,:)),'Color',cmap(2,:));

% average estimate of spectra after truncating burn-in values
simSpec_trunc = simSpec(burn_in:end,:);
simSpec_mean = mean(abs(simSpec_trunc));
plot(f, abs(simSpec_mean),'Color',cmap(1,:));


% % initial estimate of crater geometry
% shape0 = geomFunction(geomParams, geomStyle);
% figure(4); clf;
% depth0 = shape0(1,1);
% plot(shape0(:,2), shape0(:,1),'Color','k');
% set(gca,'YDir','Reverse'); hold on;
% plot(-shape0(:,2), shape0(:,1),'Color','k');
% plot([-shape0(1,2) shape0(1,2)],[depth0 depth0],'Color','k')
% 
% % final estimate of crater geometry


% figure(4); clf;
% plot(data_freq(:,1), data_freq(:,2),'k');
% xlim([0 5]);
% hold 
% 
% plot(f, abs(simSpec(1,:)));
% plot(f, abs(simSpec(end,:)));