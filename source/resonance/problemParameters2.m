function M = problemParameters2(gamma,R,TA,TC)
%
% Store properties of atmosphere and crater in structure M
% load temperatures in C

M.gamma = gamma; % ratio of heat capacities
M.r = 1000; % distance from vent to observer to calculate excess pressure at
M.R = R; % specific gas constant for dry air

% ambient pressure 
% this is assumed to be constant and set equal to atmospheric pressure
M.p = 101325; % atmospheric pressure [Pa]
M.pC = M.p;

% atmospheric properties
M.TA = TA+273.15; % temperature of atmosphere [K]
M.rhoA = M.p/(M.R*M.TA);  % density [kg/m^3]
M.cA = sqrt(M.gamma*M.R*M.TA); % speed of sound [m/s]

% properties inside crater
M.TC = TC+273.15; % temperature inside crater
M.rhoC = M.p/(M.R*M.TC);  % density [kg/m^3]
M.cC = sqrt(M.gamma*M.R*M.TC); % speed of sound [m/s]


