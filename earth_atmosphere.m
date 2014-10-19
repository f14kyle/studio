function [T P rho nu V_sonic] = earth_atmosphere(z)

%{
Earth Atmosphere

Computes the International Standard Atmosphere based on a function of 
geopotential altitude.

INPUTS
    z        - Geopotential Altitude [m]

OUTPUTS
    T        - Temperature [K]
    P        - Pressure [Pa]
    rho      - Density [kg/m^3]
    nu       - Kinematic Viscosity [m^2/s]
    V_sonic  - Speed of sound [m/s]

REFERENCE
http://aerosol.ucsd.edu/classes/sio217a/SIO217ALecture02aLMR.pdf
http://web.me.com/gyatt/atmosculator/The%20Standard%20Atmosphere.html

CONTACT
Kyle Tsai
f14kyle.work@gmail.com

UPDATED 
10/31/2011
%}

%% CONSTANTS
%  retrieved from
%  http://aerosol.ucsd.edu/classes/sio217a/SIO217ALecture02aLMR.pdf
g       = 9.80665;                  %  Gravitational Acceleration [m/s^2]
Ru      = 8.31432;                  %  Universal Gas Constant [J/(mol*K)]
MM_air  = 0.0289644;                %  Molar Mass of Dry Air [kg/mol]
R       = Ru/MM_air;                %  Gas Constant for Air 28 [J/kg-K]
gamma   = 1.4;                      %  Cp/Cv

%  VISCOSITY (EMPIRICAL CONSTANTS)
S = 110.4;                          %  Sutherland Constant
b = 0.000001458;

%  GEOPOTENTIAL HEIGHT LAYERS [km]
h0 = 0.000;                         %  Troposphere
h1 = 11000;                         %  Tropopause
h2 = 20000;                         %  Stratosphere
h3 = 32000;                         %  Stratosphere
h4 = 47000;                         %  Stratosphere
h5 = 51000;                         %  Mesosphere
h6 = 71000;                         %  MesosphereH
h7 = 84852;                         %  Mesopause

%  LAPSE RATES [degs K/km]
L0 = -6.5/1000;
L1 = 0.0/1000;
L2 = 1.0/1000;
L3 = 2.8/1000;
L4 = 0.0/1000;
L5 = -2.8/1000;
L6 = -2.0/1000;

%  BASE TEMPERATURE [degs K]
T0 = 15 + 273.15;
T1 = -56.5 + 273.15;
T2 = -56.5 + 273.15;
T3 = -44.5 + 273.15;
T4 = -2.5 + 273.15;
T5 = -2.5 + 273.15;
T6 = -58.5 + 273.15;
T7 = -86.2 + 273.15;

%  BASE ATMOSPHERIC PRESSURE [Pa]
P0 = 101325;
P1 = 22632;
P2 = 5474.9;
P3 = 868.02;
P4 = 110.91;
P5 = 66.939;
P6 = 3.9564;
P7 = 0.3734;

%% CALCULATIONS
%  retrieved from
%  http://web.me.com/gyatt/atmosculator/The%20Standard%20Atmosphere.html

%  GEOPOTENTIAL ALTITUDE [km]
h = z;                              

%   LAYER 0-1: TROPOSPHERE 0 - 11 KILOMETERS
%   Temperature decreases with altitude
if h <= h1
    TR0 = T0/T0;
    TR = TR0 + h*L0/T0;
    PR = TR^(g/(-L0*R));
end

%   LAYER 1-2: ISOTHERMAL
if h > h1 && h <= h2
    TR1 = T1/T0;
    PR1 = P1/P0;
    TR = TR1;
    PR = PR1*exp(-(h-h1)*g/(R*T1));
end

%   LAYER 2-3: INVERSION
if h > h2 && h <= h3
    TR2 = T2/T0;
    PR2 = P2/P0;
    TR = TR2 + (h - h2)*L2/T0;
    PR = PR2*(TR/TR2)^(g/(-L2*R));
end

%   LAYER 3-4: INVERSION
if h > h3 && h <= h4
    TR3 = T3/T0;
    PR3 = P3/P0;
    TR = TR3 + (h - h3)*L3/T0;
    PR = PR3*(TR/TR3)^(g/(-L3*R));
end

%   LAYER 4-5: ISOTHERMAL
if h > h4 && h <= h5
    TR4 = T4/T0;
    PR4 = P4/P0;
    TR = TR4;
    PR = PR4*exp(-(h-h4)*g/(R*T4));
end

%   LAYER 5-6
if h > h5 && h <= h6
    TR5 = T5/T0;
    PR5 = P5/P0;
    TR = TR5 + (h - h5)*L5/T0;
    PR = PR5*(TR/TR5)^(g/(-L5*R));
end

%   LAYER 6-7
if h > h6 && h <= h7
    TR6 = T6/T0;
    PR6 = P6/P0;
    TR = TR6 + (h - h6)*L6/T0;
    PR = PR6*(TR/TR6)^(g/(-L6*R));
end

T           = TR*T0;                        %  Local Temperature [K]
P           = PR*P0;                        %  Local Pressure [Pa]
rho         = P/(R*T);                      %  Local Density [kg/m^3]
mu          = (b*T^1.5/(S+T));              %  Dynamic Viscosity [kg/m-s]
nu          = mu/rho;                       %  Kinematic Viscosity [m^2/s]
V_sonic     = (gamma*R*T)^0.5;              %  Speed of Sound [m/s]

end