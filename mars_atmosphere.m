% aero_getMarsAtmo.m

% DESCRIPTION

% Uses data taken from Sample Mars-GRAM 2005 Auxiliary Profile (Average of 
% 17 TES Limb Profiles at Phoenix Landing conditions.  Actual values 
% retrieved from Table 1 of Atmospheric Models for Mars Aerocapture by 
% C.G. Justus, Aleta Duvall, and Vernon W. Keller
% http://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20050207555.pdf

% Ideal gas law is used to compute densities of helium and hydrogen on
% Mars.  

% http://pds-atmospheres.nmsu.edu/education_and_outreach/encyclopedia/gas_constant.htm

% INPUTS
% h        - altitude [km]

% OUTPUTS
% rho_amb      - density of ambient air [kg/m^3]
% P            - pressure [Pa]
% T            - temperature [K]
% U_sonic      - speed of sound [m/s]
% rho_H2       - hydrogen density [kg/m^3]
% rho_He       - helium density [kg/m^3]

% WRITTEN BY
% Kyle Tsai
% ktsai6@illinois.edu

% LAST MODIFIED
% 3/09/2014

%{
MARS FACTS
Surface pressure:  6.36 mb at mean radius (variable from 4.0 to 8.7 mb depending on season)  
                   [6.9 mb to 9 mb (Viking 1 Lander site)]
Surface density: ~0.020 kg/m3
Scale height:  11.1 km
Total mass of atmosphere: ~2.5 x 1016 kg
Average temperature:  ~210 K (-63 C)
Diurnal temperature range: 184 K to 242 K (-89 to -31 C) (Viking 1 Lander site)
Wind speeds:  2-7 m/s (summer), 5-10 m/s (fall), 17-30 m/s (dust storm) (Viking Lander sites)
Mean molecular weight: 43.34 g/mole 
Atmospheric composition (by volume): 
    Major      : Carbon Dioxide (CO2) - 95.32% ; Nitrogen (N2) - 2.7%
                 Argon (Ar) - 1.6%; Oxygen (O2) - 0.13%; Carbon Monoxide (CO) - 0.08% 
    Minor (ppm): Water (H2O) - 210; Nitrogen Oxide (NO) - 100; Neon (Ne) - 2.5;
                 Hydrogen-Deuterium-Oxygen (HDO) - 0.85; Krypton (Kr) - 0.3; 
		 Xenon (Xe) - 0.08

http://nssdc.gsfc.nasa.gov/planetary/factsheet/marsfact.html
%}


function [P T nu U_sonic rho_amb rho_He rho_H2] = getMarsAtmo(h)

% GIVEN
% Vector of heights [km]
h_vector = [2.172 4.825 7.408 9.921 12.365 14.744 17.062 19.326 21.54 23.709 25.838 27.931 29.992 32.029 34.048 36.057 38.061 40.051 42.039 44.019 45.987 47.948 49.888 51.825 53.72 55.646 57.546];
% Vector of temperatures [K]
T_vector = [209.39 203.59 197.79 192.08 186.63 181.46 176.77 172.53 168.64 165.22 162.06 159.2 156.83 155.03 153.86 153.15 152.5 151.68 150.77 149.83 148.78 147.71 146.72 145.76 144.75 143.68 142.68];
% Vector of pressures [Pa]
P_vector = [4.7513E02 3.7003E02 2.8813E02 2.2443E02 1.7483E02 1.3613E02 1.0603E02 8.2553E01 6.4293E01 5.0073E01 3.9003E01 3.0373E01 2.3653E01 1.8423E01 1.4353E01 1.1173E01 8.7003E00 6.7803E00 5.2803E00 4.1103E00 3.2003E00 2.4903E00 1.9403E00 1.5103E00 1.1803E00 9.1703E-01 7.1403E-01];
% Vector of densities [kg/m^3]
rho_vector = [1.1873E-02 9.5093E-03 7.6223E-03 6.1133E-03 4.9003E-03 3.9253E-03 3.1383E-03 2.5043E-03 1.9953E-03 1.5863E-03 1.2593E-03 9.9823E-04 7.8913E-04 6.2183E-04 4.8813E-04 3.8223E-04 2.9913E-04 2.3443E-04 1.8373E-04 1.4393E-04 1.1293E-04 8.8513E-05 6.9463E-05 5.4463E-05 4.2893E-05 3.3603E-05 2.6363E-05];

% INTERPOLATATION
P       = interp1(h_vector,P_vector,h,'linear','extrap');
T       = interp1(h_vector,T_vector,h,'linear','extrap');
rho_amb = interp1(h_vector,rho_vector,h,'linear','extrap');

% CALCULATION
% Calculate densities of helium and hydrogen using the ideal gas law.
MM_He = 4.0026022/1000;         % molar mass of helium [kg/mol]
MM_H2 = 2.015894/1000;          % molar mass of hydrogen [kg/mol]
MM_amb = 44.01/1000;            % molar mass of ambient Mars atmo [kg/mol]
R_u = 8.3143;                   % univeral gas constant [J/mol-K]
R = R_u/MM_amb;                 % specific gas constant [J/kg-K]
rho_He = P * MM_He / (R * T);   % density of helium [kg/m^3]
rho_H2 = P * MM_H2 / (R * T);   % density of hydrogen [kg/m^3]

% Calculate the speed of sound
gamma = 1.2941;                 % ratio of Cp to Cv
U_sonic = (gamma * R * T)^0.5;      % speed of sound [m/s] 

% Calculate the dynamic viscosity using Sutherland's formula
C = 240;                        % Sutherland's constant for CO2
T_0 = 293.15;                   % C02 reference temperature [K]
mu_0 = 14.8E-6;                 % CS02 reference viscosity
mu = mu_0*(T_0 + C)/(T + C)*(T/T_0)^1.5;    % dynamic viscosity (kg/m-s)
nu = mu/rho_amb;                % kinematic viscosity [m2/s]