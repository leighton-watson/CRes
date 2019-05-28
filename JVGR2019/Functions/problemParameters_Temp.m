% OLD: M = problemParameters()
% NEW: no longer a function
function M = problemParameters_Temp(T)
%
% Store properties of atmosphere and crater in structure M

M.gamma = 1.4; % ratio of heat capacities
M.r = 1000; % distance from vent to observer to calculate excess pressure at
M.R = 287; % specific gas constant for dry air

%%% updated version (October 12, 2018)

% temperature
M.TA = 0; % temperature of atmosphere [C]
M.TC = T; % temperature inside crater [C]

% speed of sound in atmosphere
M.cA = sqrt(M.gamma*M.R*(M.TA+273.15)); % [m/s]

% density of atmospheric air
M.rhoA = 1; % [kg/m^3]

% pressure in atmosphere
M.pA = M.cA^2.*M.rhoA/M.gamma; % [Pa]

% background pressure in crater is equal to pressure in atmosphere
M.pC = M.pA; % [Pa]

% speed of sound in crater
M.cC = sqrt(M.gamma*M.R*(M.TC+273.15));

% density of air inside crater
M.rhoC = M.pC.*M.gamma./(M.cC.^2);

% % atmospheric properties
% M.TA = 0; % temperature of atmosphere
% M.rhoA = 1; % density [kg/m^3]
% M.cA = sqrt(M.gamma*M.R*(M.TA+273.15)); % speed of sound [m/s]
% M.pA = M.cA^2*M.rhoA/M.gamma; % pressure [Pa]
% 
% % properties inside crater
% M.TC = T+273.15; % temperature inside crater
% M.rhoC = M.pA./(M.R*M.TC); % density inside crater
% %M.rhoC = 1; % density [kg/m^3]
% %M.cC = sqrt(M.gamma*M.R*(M.TC+273.15)); % speed of sound [m/s]
% M.cC = sqrt(M.gamma*M.R*(M.TC));
% M.pC = M.pA*ones(size(T)); %M.cC.^2.*M.rhoC./M.gamma; % pressure [Pa]
% 
