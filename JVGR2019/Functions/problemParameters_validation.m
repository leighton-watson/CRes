function M = problemParameters_validation(gamma, r)
% M = problemParameters(gamma)
%
% Store properties of atmosphere and crater in structure M

M.gamma = gamma; % ratio of heat capacities
M.r = r; % distance from vent to observer to calculate excess pressure at
M.R = 287.06; % specific gas constant for dry air

% atmospheric properties
M.TA = 0; % temperature of atmosphere
M.rhoA = 1; % density [kg/m^3]
M.cA = sqrt(M.gamma*M.R*(M.TA+273.15)); % speed of sound [m/s]

% properties inside crater
M.TC = 0; % temperature inside crater
M.rhoC = 1; % density [kg/m^3]
M.cC = sqrt(M.gamma*M.R*(M.TC+273.15)); % speed of sound [m/s]
M.pC = M.cC^2*M.rhoC/M.gamma; % pressure [Pa]

