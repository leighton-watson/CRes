%% DISP RESULTS %%
%
% Display inversions results

clear all;
clc;
path(pathdef);
cmap = get(gca,'ColorOrder');

addpath Data/Fig12_Inversion/
addpath Functions/
addpath Inversion

figure(1); clf;

%% DATA GEOMETRY %%


%%% Structure from motion %%%
SfM = load('SfM_Geom.csv');

figure(1); 
subplot(2,3,1);
hold on;
box on;

SfM_avg = mean(SfM(:,2:5),2);

plot(SfM_avg, SfM(:,1),'Color','k');
plot(-SfM_avg, SfM(:,1),'Color','k');

% for i = 1:4
%     plot(SfM(:,i+1),SfM(:,1),'Color','k','LineStyle','-');
%     plot(-SfM(:,i+1),SfM(:,1),'Color','k','LineStyle','-');
% end

set(gca,'YDir','Reverse');
xlabel('Radius (m)');
ylabel('Depth (m)');
xlim([-110 110])
ylim([0 160])

subplot(2,3,4);
hold on;
box on;

plot(SfM_avg, SfM(:,1),'Color','k');
plot(-SfM_avg, SfM(:,1),'Color','k');

set(gca,'YDir','Reverse');
xlabel('Radius (m)');
ylabel('Depth (m)');
xlim([-110 110])
ylim([0 160])


%% EARLY TIME %%
load output_early

% plot spectra data
subplot(2,3,[2 3]);
plot(data_freq(:,1), abs(data_freq(:,2))./max(abs(data_freq(:,2))),'k');
xlabel('Frequency (Hz)');
xlim([0 4]);
hold on;
ylim([0 1.2]);

% resonant properties of data
data.f = data_freq(:,1);
data.P = data_freq(:,2);
[f0d1,Qd1] = resPeakProps(data,'data')

% plot spectra
%plot(f, abs(simSpec_mean),'Color',cmap(1,:));

% % initial estimate of spectra
% plot(f, abs(simSpec(1,:)),'Color',cmap(3,:));

flow = 0.2;
fhigh = 0.5;

% create filter
for i = 1:length(f)
    if f(i) < flow
        spec_filt(i) = 0;
    elseif f(i) < fhigh
        spec_filt(i) = (0.5*(sin(pi*(f(i)-flow-(fhigh-flow)/2)/(fhigh-flow))+1)).^3;
    else
        spec_filt(i) = 1;
    end
end


%figure(3); clf;
%plot(f, spec_filt);

% plot filtered spectra
figure(1); subplot(2,3,[2 3]);
plot(f, abs(simSpec_mean).*spec_filt,'Color',cmap(1,:));

% resonant properties of average spectra
sim_mean.f = f;
sim_mean.P = simSpec_mean;
[f0s1mean, Qs1mean] = resPeakProps(sim_mean, 'sim')

plot(f, abs(simSpec(1,:)).*spec_filt./max(abs(simSpec(1,:)).*spec_filt),'Color',cmap(3,:));

% resonant properties of input
sim_input.f = f;
sim_input.P = simSpec(1,:);
[f0s1input, Qs1input] = resPeakProps(sim_input,'sim')


% standard deviation of spectra
simSpec_std = std(abs(simSpec_trunc));
plot(f, (abs(simSpec_mean)+abs(simSpec_std)).*spec_filt,'Color',cmap(1,:),'LineStyle',':');
plot(f, (abs(simSpec_mean)-abs(simSpec_std)).*spec_filt,'Color',cmap(1,:),'LineStyle',':');

% plot crater geometry
subplot(2,3,1);
plot(shapeF(:,2), shapeF(:,1),'Color',cmap(1,:));
plot(-shapeF(:,2), shapeF(:,1),'Color',cmap(1,:));
plot([-shapeF(1,2) shapeF(1,2)],[depthF depthF],'Color',cmap(1,:));

% standard deviation of geometry
x_std = std(x_trunc);
geomParams_mean_stdp = (x_mean+x_std).*geomParams;
shape_stdp = geomFunction(geomParams_mean_stdp,geomStyle);
depth_stdp = shape_stdp(1,1);
plot(shape_stdp(:,2), shape_stdp(:,1),'Color',cmap(1,:),'LineStyle',':');
plot(-shape_stdp(:,2), shape_stdp(:,1),'Color',cmap(1,:),'LineStyle',':');
plot([-shape_stdp(1,2) shape_stdp(1,2)],[depth_stdp depth_stdp],...
    'Color',cmap(1,:),'LineStyle',':');

geomParams_mean_stdm = (x_mean-x_std).*geomParams;
shape_stdm = geomFunction(geomParams_mean_stdm,geomStyle);
depth_stdm = shape_stdm(1,1);
plot(shape_stdm(:,2), shape_stdm(:,1),'Color',cmap(1,:),'LineStyle',':');
plot(-shape_stdm(:,2), shape_stdm(:,1),'Color',cmap(1,:),'LineStyle',':');
plot([-shape_stdm(1,2) shape_stdm(1,2)],[depth_stdm depth_stdm],...
    'Color',cmap(1,:),'LineStyle',':');

% plot initial estimate of crater geometry
plot(shapeI(:,2), shapeI(:,1),'Color',cmap(3,:));
plot(-shapeI(:,2), shapeI(:,1),'Color',cmap(3,:));
plot([-shapeI(1,2) shapeI(1,2)],[depthI depthI],'Color',cmap(3,:));

%%% compute resonant properties %%%





%% LATE TIME %%
load output_late

% plot spectra data
figure(1);
subplot(2,3,[5 6]);
plot(data_freq(:,1), abs(data_freq(:,2))./max(abs(data_freq(:,2))),'k');
xlabel('Frequency (Hz)');
xlim([0 4]);
hold on;
ylim([0 1.2])

% resonant properties of data
data.f = data_freq(:,1);
data.P = data_freq(:,2);
[f0d2,Qd2] = resPeakProps(data,'data')

% plot spectra
%plot(f, abs(simSpec_mean),'Color',cmap(2,:));

% % plot filtered spectra
% plot(f, abs(simSpec_mean).*spec_filt,'Color',cmap(4,:));
% 
% % initial estimate of spectra
% plot(f, abs(simSpec(1,:)),'Color',cmap(3,:));

% plot filtered spectra
figure(1); subplot(2,3,[5 6]);
plot(f, abs(simSpec_mean).*spec_filt,'Color',cmap(2,:));

% resonant properties of average spectra
sim_mean.f = f;
sim_mean.P = simSpec_mean;
[f0s2mean, Qs2mean] = resPeakProps(sim_mean, 'sim')

plot(f, abs(simSpec(1,:)).*spec_filt./max(abs(simSpec(1,:)).*spec_filt),'Color',cmap(3,:));

% resonant properties of input
sim_input.f = f;
sim_input.P = simSpec(1,:);
[f0s2input, Qs2input] = resPeakProps(sim_input,'sim')

% standard deviation of spectra
simSpec_std = std(abs(simSpec_trunc));
plot(f, (abs(simSpec_mean)+abs(simSpec_std)).*spec_filt,'Color',cmap(2,:),'LineStyle',':');
plot(f, (abs(simSpec_mean)-abs(simSpec_std)).*spec_filt,'Color',cmap(2,:),'LineStyle',':');


% plot crater geometry
subplot(2,3,4);
plot(shapeF(:,2), shapeF(:,1),'Color',cmap(2,:));
plot(-shapeF(:,2), shapeF(:,1),'Color',cmap(2,:));
plot([-shapeF(1,2) shapeF(1,2)],[depthF depthF],'Color',cmap(2,:));

% standard deviation of geometry
x_std = std(x_trunc);
geomParams_mean_stdp = (x_mean+x_std).*geomParams;
shape_stdp = geomFunction(geomParams_mean_stdp,geomStyle);
depth_stdp = shape_stdp(1,1);
plot(shape_stdp(:,2), shape_stdp(:,1),'Color',cmap(2,:),'LineStyle',':');
plot(-shape_stdp(:,2), shape_stdp(:,1),'Color',cmap(2,:),'LineStyle',':');
plot([-shape_stdp(1,2) shape_stdp(1,2)],[depth_stdp depth_stdp],...
    'Color',cmap(2,:),'LineStyle',':');

geomParams_mean_stdm = (x_mean-x_std).*geomParams;
shape_stdm = geomFunction(geomParams_mean_stdm,geomStyle);
depth_stdm = shape_stdm(1,1);
plot(shape_stdm(:,2), shape_stdm(:,1),'Color',cmap(2,:),'LineStyle',':');
plot(-shape_stdm(:,2), shape_stdm(:,1),'Color',cmap(2,:),'LineStyle',':');
plot([-shape_stdm(1,2) shape_stdm(1,2)],[depth_stdm depth_stdm],...
    'Color',cmap(2,:),'LineStyle',':');

% plot initial estimate of crater geometry
plot(shapeI(:,2), shapeI(:,1),'Color',cmap(3,:));
plot(-shapeI(:,2), shapeI(:,1),'Color',cmap(3,:));
plot([-shapeI(1,2) shapeI(1,2)],[depthI depthI],'Color',cmap(3,:));


