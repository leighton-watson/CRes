function [x, f0, Q, misfit, simSpec, f, count] = mcmc(nIter, dx, geomParams, geomFlag, geomStyle, srcParams, srcFlag, srcStyle, lowerBnds, upperBnds, discrParams, misfitParams, data);
% [x f0 Q misfit simSpec f count] = mcmc(nIter dx geomParams geomFlag srcParams srcFlag lowerBnds upperBnds discParams data)
% 
% Inputs:
% nIter = number of iterations
% dx = step size 
% geomParams = parameters of crater geometry
% geomFlag = Boolean that specifies if geometry parameters are inverted for or are kept constant
% srcParams = parameters of source function
% srcFlag = Boolean that specifies if source parameters are inverted for or are kept constant
% lowerBnds = lower bounds of parameters that are inverted for
% upperBnds = upper bounds of parameters that are inverted for
% discrParms = parameters of discretization 
% data = observed amplitude spectrum
%
% Outputs:
% x = inverted parameters at each successful iteration 
% f0 = simulated resonant frequency at each successful iteration
% Q = simulated quality factor at each successful iteration
% misfit = misfit at each successful iteration
% simSpec = simulated amplitude spectrum at each successful iteration
% f = frequency vector that corresponds to simSpec
% count = number of successful iterations. count/nIter = acceptance rate

%% Discretization Information and Resonance1d Params %%

N = discrParams(2); % number of grid points
Nf = discrParams(3); % number of frequency samples
freq = [0 discrParams(4)]; % frequency vector (Hz)
order = 4; % order of numerical scheme (4, 6 or 8)
style = 'baffled piston'; % acoustic radiation model ('monopole' or ' baffled piston')
M = problemParameters(); % problem parameters required for resonance1d

%% Load Variables %%

geomLength = length(geomParams); % length of geomParams vector
srcLength = length(srcParams); % length of srcParams vector

if ~geomFlag && ~srcFlag
    error('Both geomFlag and srcFlag are set equal to zero. There are no parameters to invert for.')
end 

vars = []; % variables that are inverted for
if geomFlag % include geometry parameters in inversion
    vars = [vars geomParams];
end
if srcFlag % include source parameters in inversion
    vars = [vars srcParams];
end

x0 = ones(size(vars)); % vector of normalized inversion variables

%% Initialize MCMC %%

count = 0; % number of successful iterations. count/nIter = acceptance rate
alpha = 0.9; % steps that increase misfit are accepted if a randomly generated number exceeds this value
beta = misfitParams(1); % parameter in the objective function
delta = misfitParams(2); % parameter in the objective function
gamma = misfitParams(3); % parameter in the objective function

%% Initial Inversion

x(1,:) = x0; % initial estimate of inversion variables

% compute geometry
if geomFlag 
    shape = geomFunction(x(1,1:geomLength).*vars(1:geomLength), geomStyle);
else
    shape = geomFunction(geomParams, geomStyle);
end
depth = shape(1,1);

% compute source
if srcFlag
    if geomFlag
        [S, ~, ~, ~] = sourceFunction(1, x(1,geomLength+1:geomLength+srcLength).*vars(geomLength+1:geomLength+srcLength), srcStyle, discrParams);
    else
        [S, ~, ~, ~] = sourceFunction(1, x(1,1:srcLength).*vars(1:srcLength), srcStyle, discrParams);
    end        
else
    [S, ~, ~, ~] = sourceFunction(1, srcParams, srcStyle, discrParams);
end

% simulate spectra
res = resonance1d(shape, depth, freq, Nf, style, order, M); % compute transfer function
sim.f = res.f; % frequency vector
sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1)); % convolve transfer function and source function and normalize amplitudes
 
% data
dataSpectrum.P = pchip(data(:,1), data(:,2), sim.f); % interpolate data on same regular frequency vector as simulations
dataSpectrum.P = dataSpectrum.P./max(dataSpectrum.P); % normalize amplitude
dataSpectrum.f = sim.f;

% compute resonant frequency and quality factor for simulation and data
[data_f0, data_Q] = resPeakPropsInversion(dataSpectrum,'data');
[sim_f0, sim_Q] = resPeakPropsInversion(sim,'sim');

% compute misfit
misfit(1) = delta*abs(data_f0-sim_f0)/data_f0 + gamma*abs(data_Q - beta*sim_Q)/data_Q;

% save outputs
simSpec(1,:) = sim.P;
f0(1) = sim_f0;
Q(1) = sim_Q;
f = sim.f;

%% Iterate %%

for i = 2:nIter
   
    if mod(i,10) == 0 % print iteration number
        disp(strcat('Iteration Number:',num2str(i),'/',num2str(nIter)));
    end
    
    satBnds = false; % Boolean that checks if upper and lower bounds are satisfied
    
    % save values from previous iteration
    x_latest = x(end,:);
    misfit_latest = misfit(end,:);
    
    % check bounds
    while satBnds == false
        r = -1 + 2*rand(size(x0)); % generate vector of random numbers that is the same length as the number inversion variables
        x_potential = x_latest + dx*r; % take random step in inversion variables
        rBnds = x_potential.*vars; % compute physical value of inversion variables (non-normalized)
        
        if strcmp(geomStyle,'intp')
            if sum(rBnds>lowerBnds)==length(vars) && sum(rBnds<upperBnds)==length(vars) && ... 
                all(diff(x_potential(2:geomLength).*geomParams(2:end))<0)
                satBnds = true;
            end
            
        else
            if sum(rBnds>lowerBnds)==length(vars) && sum(rBnds<upperBnds)==length(vars) % check if inverted variables satisfy bounds
                satBnds = true;
            end
        
        end
        
    end
    
   
    % compute geometry
    if geomFlag % if geometry variables are updated compute new geometry
        shape = geomFunction(x_potential(1:geomLength).*vars(1:geomLength), geomStyle);
        depth = shape(1,1);
    end
    
    % compute source
    if srcFlag
        if geomFlag
            [S, ~, ~, ~] = sourceFunction(1, x_latest(geomLength+1:geomLength+srcLength).*vars(geomLength+1:geomLength+srcLength), srcStyle, discrParams);
        else
            [S, ~, ~, ~] = sourceFunction(1, x_latest(1:srcLength).*vars(1:srcLength), srcStyle, discrParams);
        end
    end
            
    
    % simulate spectra
    res = resonance1d(shape, depth, freq, Nf, style, order, M); % compute transfer function
    sim.f = res.f; % frequency vector
    sim.P = (res.P(1:N/2+1).*S(1:N/2+1))./max(res.P(1:N/2+1).*S(1:N/2+1)); % convolve transfer function and source function and normalize amplitudes
 
    % compute resonant frequency and quality factor for simulation and data
    [sim_f0, sim_Q] = resPeakPropsInversion(sim,'sim');

    % compute misfit
    misfit_potential = delta*abs(data_f0-sim_f0)/data_f0 + gamma*abs(data_Q - beta*sim_Q)/data_Q;     
    
    if misfit_potential < misfit_latest % misfit has decreased
        x = [x; x_potential];
        misfit = [misfit; misfit_potential];
        simSpec = [simSpec; sim.P];
        f0 = [f0; sim_f0];
        Q = [Q; sim_Q];
        count = count + 1;
        
    elseif rand(1) > alpha % stochastic process that accepts some steps that increase misfit
        x = [x; x_potential];
        misfit = [misfit; misfit_potential];
        simSpec = [simSpec; sim.P];
        f0 = [f0; sim_f0];
        Q = [Q; sim_Q];
        count = count + 1;
        
    end
    
end





