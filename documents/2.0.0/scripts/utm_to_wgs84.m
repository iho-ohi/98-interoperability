function [lat, lon] = utm_to_wgs84(E, N, Nz, is_north)
% UTM_TO_WGS84 Convert UTM coordinates to WGS84 latitude/longitude
%
% Inputs:
%   E        - Easting (meters)
%   N        - Northing (meters)
%   Nz       - UTM zone number (1-60)
%   is_north - true for North zones, false for South zones
%
% Outputs:
%   lat      - Latitude in degrees
%   lon      - Longitude in degrees
%
% Based on Appendix H-1.1 of S-100 ECDIS standard

    % Constants for all UTM zones
    FE = 500000;  % False easting (m)
    k0 = 0.9996;  % Scale at natural origin
    a = 6378137.000;  % Semi-major axis (m)
    inv_f = 298.257223563;  % Inverse flattening

    % Calculate derived constants
    f = 1 / inv_f;  % = 3.35281066475e-3
    b = a * (1 - f);  % Semi-minor axis = 6356752.31425

    % False northing
    if is_north
        FN = 0;
    else
        FN = 10000000;
    end

    % Calculate central meridian λ0
    if Nz >= 1 && Nz <= 30
        lambda0 = 183 - (6 * Nz);  % degrees West (negative)
        lambda0 = -lambda0;  % Convert to East (positive East)
    elseif Nz >= 31 && Nz <= 60
        lambda0 = (6 * (Nz - 30)) - 3;  % degrees East
    else
        error('Zone number must be between 1 and 60');
    end

    % Calculate projection constants
    e = sqrt(2*f - f^2);  % = 8.18191908e-2
    n = f / (2 - f);  % = 1.67922039e-3
    B = (a / (1+n)) * (1 + n^2/4 + n^4/64);  % = 6367449.145823415

    h1_prime = n/2 - (2/3)*n^2 + (37/96)*n^3 + (1/360)*n^4;  % = 8.37732164e-4
    h2_prime = (1/48)*n^2 - (1/15)*n^3 + (437/1440)*n^4;  % = 5.84321837e-8
    h3_prime = (17/480)*n^3 - (37/840)*n^4;  % = 1.67348888e-10
    h4_prime = (4397/161280)*n^4;  % = 2.16773776e-13

    % Step 5: Calculate η' and ξ'
    eta_prime = (E - FE) / (B * k0);
    xi_prime = (N - FN) / (B * k0);

    % Calculate ξ'_i and η'_i terms
    xi1_prime = h1_prime * sin(2*xi_prime) * cosh(2*eta_prime);
    eta1_prime = h1_prime * cos(2*xi_prime) * sinh(2*eta_prime);

    xi2_prime = h2_prime * sin(4*xi_prime) * cosh(4*eta_prime);
    eta2_prime = h2_prime * cos(4*xi_prime) * sinh(4*eta_prime);

    xi3_prime = h3_prime * sin(6*xi_prime) * cosh(6*eta_prime);
    eta3_prime = h3_prime * cos(6*xi_prime) * sinh(6*eta_prime);

    xi4_prime = h4_prime * sin(8*xi_prime) * cosh(8*eta_prime);
    eta4_prime = h4_prime * cos(8*xi_prime) * sinh(8*eta_prime);

    xi0_prime = xi_prime - (xi1_prime + xi2_prime + xi3_prime + xi4_prime);
    eta0_prime = eta_prime - (eta1_prime + eta2_prime + eta3_prime + eta4_prime);

    beta_prime = asin(sin(xi0_prime) / cosh(eta0_prime));

    Q_prime = asinh(tan(beta_prime));

    % Iterate to find Q"
    Q_double_prime = Q_prime + e * atanh(e * tanh(Q_prime));

    tolerance = 1e-8;
    max_iterations = 10;
    for iter = 1:max_iterations
        Q_old = Q_double_prime;
        Q_double_prime = Q_prime + e * atanh(e * tanh(Q_double_prime));

        if abs(Q_double_prime - Q_old) < tolerance
            break;
        end
    end

    % Step 6: Calculate final latitude and longitude
    lat = atan(sinh(Q_double_prime)) * (180/pi);
    lon = lambda0 + asin(tanh(eta0_prime) / cos(beta_prime)) * (180/pi);
end
