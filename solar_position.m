function [solar_hour declination elevation zenith azimuth] = ...
    solar_position(day,hour,minute,latitude,longitude,dGMT)

%{
Solar Position Calculator

This script calculates the sun's position as described by the horizontal
and equatorial coordinate systems as a function of time and location.

INPUTS
    day - [1-366]
    hour         - [1-24]
    minute       - [1-59]
    latitude     - [degs]
    longitude    - [degs]
    dGMT         - Difference of local time from GMT [hrs]

OUTPUTS
    solar_hour   - Solar Hour Angle [degs]
    declination  - Declination Angle [degs]
    elevation    - Elevation Angle [degs]
    azimuth      - Azimuth Angle [degs]
    zenith       - Zenith Angle [degs]
    
REFERENCE
    http://pvcdrom.pveducation.org/SUNLIGHT/SOLART.HTM

CONTACT
    Kyle Tsai
    f14kyle.work@gmail.com

UPDATED
    10/31/2011
    %}
    
    %  LOCAL STANDARD TIME MERIDIAN (LSTM) [degs]
    %  LSTM is a reference meridian used for a particular time zone.
    LSTM = 15 * dGMT;
    
    %  EQUATION OF TIME (EoT)[degs]
    %  EoT corrects for the eccentricity of the orbit and the axial tilt.
    B = (360/365)*(day - 81);
    EoT = 9.87*sind(2*B) - 7.53*cosd(B) - 1.5*sind(B);
    
    %  TIME CORRECTION FACTOR (TC) [min]
    %  TC accounts for the variation of the Local Solar Time (LST) within a
    %  given time zone due to the longitude variation within the time zone.
    TC = 4*(longitude - LSTM) + EoT;
    
    %  LOCAL SOLAR TIME (LST) [hrs]
    LST = hour + minute/60 + TC/60;
    
    %  SOLAR HOUR ANGLE (SHA)[degs]
    %  converts LST into number of degrees sun moves across the sky
    solar_hour = 15*(LST - 12);
    
    %  DECLINATION ANGLE [degs]
    declination = 23.45 * sin((day - 80) * (2*pi/365.25));
    
    %  ELEVATION ANGLE [degs]
    elevation = asind(sind(declination).*sind(latitude) ...
        + cosd(declination).*cosd(solar_hour).*cosd(latitude));
    
    %  ZENITH ANGLE [degs]
    zenith = 90 - elevation;
    
    %  AZIMUTH ANGLE [degs]
    if LST < 12
        azimuth = acosd((sind(declination)*cosd(latitude) - ...
            cosd(declination)*sind(latitude)*cosd(solar_hour))/cosd(elevation));
    else
        azimuth = 360 - acosd((sind(declination)*cosd(latitude) - ...
            cosd(declination)*sind(latitude)*cosd(solar_hour))/cosd(elevation));
    end
end