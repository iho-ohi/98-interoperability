function [lat, lon] = ups_north_to_wgs84(E, N)
% UPS_NORTH_TO_WGS84 Convert UPS North coordinates to WGS84 latitude/longitude
%
% Inputs:
%   E   - Easting (meters)
%   N   - Northing (meters)
%
% Outputs:
%   lat - Latitude in degrees
%   lon - Longitude in degrees
%
% Based on Appendix H-1.2 of S-100 ECDIS standard

    % Constants for UPS North
    FE = 2000000;  % False easting (m)
    FN = 2000000;  % False northing (m)
    k0 = 0.994;    % Scale at natural origin
    a = 6378137.000;  % Semi-major axis (m)
    inv_f = 298.257223563;  % Inverse flattening

    % Calculate derived constants
    f = 1 / inv_f;  % = 3.35281066475e-3
    e = sqrt(2*f - f^2);  % = 8.18191908e-2

    % Step 2: Calculate t'
    term1 = sqrt((E - FE)^2 + (N - FN)^2);
    term2 = ((1 + e)^(1+e)) * ((1 - e)^(1-e));
    t_prime = (term1 * term2) / (2 * a * k0);

    % Simplified constant calculation
    % t_prime = 7.91307065e-8 * sqrt((E - 2000000)^2 + (N - 2000000)^2);

    % Step 3: Calculate χ for UPS North
    chi = pi/2 - 2*atan(t_prime);

    % Step 4: Calculate φ using the series expansion
    % Coefficients for the series
    c1 = e^2/2 + 5*e^4/24 + e^6/12 + 13*e^8/360;  % = 3.35655147e-3
    c2 = 7*e^4/48 + 29*e^6/240 + 811*e^8/11520;  % = 6.57187270e-6
    c3 = 7*e^6/120 + 81*e^8/1120;  % = 1.76456433e-8
    c4 = 4279*e^8/161280;  % = 5.32847842e-11

    lat = (180/pi) * (chi + ...
                      c1 * sin(2*chi) + ...
                      c2 * sin(4*chi) + ...
                      c3 * sin(6*chi) + ...
                      c4 * sin(8*chi));

    % Step 5: Calculate λ
    if E == 2000000
        lon = 0;
    else
        lon = (180/pi) * atan2((E - 2000000), (2000000 - N));
    end
end
