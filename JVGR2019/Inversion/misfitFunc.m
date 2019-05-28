function misfit = misfitFunc(sim, data, source, misfitParams, resParams)
% misfit = misfitFunc(sim, data)
%
% Compute misfit function
%
% Inputs:
% sim = simulation output from resonance1d
% data = data of infrasound signal in frequency domain. First column =
% frequency, second column = amplitude
% params = weighting parameters for misfit function
%
% Outputs:
% misfit = value of misfit function

% specify parameter values
gamma = misfitParams(1);
N = resParams.N;

% convolve transfer function with source
transFunc = sim.P;
S = source;
sim_spectra = transFunc(1:N/2+1) .* S(1:N/2+1);

% Compute resonant frequency and quality factor of simulation
[fsim, Qsim] = resPeakProps(sim.f, sim_spectra);

% Compute resonant frequency and quality factor of data
[fdata, Qdata] = resPeakProps(data(:,1), data(:,2));
    
% calculate misfit function
misfit = abs(fdata-fsim) - gamma * abs(Qdata-Qsim);