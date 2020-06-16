function M = problemParametersInv(craterTemp,atmoTemp)
%
% Store properties of atmosphere and crater in structure M

M.gamma = 1.4; % ratio of heat capacities
M.r = 1000; % distance from vent to observer to calculate excess pressure at
M.R = 287; % specific gas constant for dry air

% atmospheric properties
M.TA = atmoTemp; % temperature of atmosphere
M.pA = 1e5; % pressure [Pa]
M.cA = sqrt(M.gamma*M.R*(M.TA+273.15)); % speed of sound [m/s]
M.rhoA = (M.gamma*M.pA)./(M.cA^2); % density [kg/m^3]

% properties inside crater
M.TC = craterTemp; % temperature inside crater
M.pC = 1e5; % pressure [Pa]
M.cC = sqrt(M.gamma*M.R*(M.TC+273.15)); % speed of sound [m/s]
M.rhoC = (M.gamma*M.pC)./(M.cC^2); % density [kg/m^3]


