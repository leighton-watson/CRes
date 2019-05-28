function M = problemParameters_GasComp(gamma,TA,TC,RA,RC,rhoA)
% M = problemParameters_GasComp(TA,TC,RA,RC,rhoA)
%
% Store properties of atmosphere and crater in structure M
%
% INPUTS:
% TA = atmospheric temperature [C]
% TC = crater temperature [C]
% RA = specific gas constant for atmosphere [J/kg/K]
% RC = specific gas constant for crater [J/kg/K]
% rhoA = atmospheric density [kg/m^3]


M.gamma = gamma; % ratio of heat capacities
M.r = 1000; % distance from vent to observer to calculate excess pressure at

M.TA = TA + 273.15; % atmospheric temperature
M.TC = TC + 273.15; % crater temperature

M.RA = RA; % specific gas constant for atmosphere
M.RC = RC; % specific gas constant for crater

M.rhoA = rhoA; % atmospheric density

% atmospheric speed of sound
M.cA = sqrt(M.gamma*M.RA*M.TA);

% crater speed of sound
M.cC = sqrt(M.gamma*M.RC*M.TC);

% atmospheric pressure
M.pA = M.cA.^2.*M.rhoA./M.gamma;

% crater pressure is the same as atmospheric pressure
M.pC = M.pA;

% calculate crater density
M.rhoC = M.pC.*M.gamma./(M.cC.^2);



