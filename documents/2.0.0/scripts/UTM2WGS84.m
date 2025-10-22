function [lat, lon] = UTM2WGS84(easting, northing)
% UTM2WGS84 Convert projected coordinates to WGS84 latitude/longitude
% This function translates the reverse() method from Proj.java
%
% Inputs:
%   easting  - Easting coordinate (meters)
%   northing - Northing coordinate (meters)
%
% Outputs:
%   lat - Latitude in degrees
%   lon - Longitude in degrees
%
% Note: This implementation uses the ellipsoid parameters from the original
% Java code (a=6377563.396, 1/f=299.3249646) which appears to be OSGB36,
% not WGS84. The projection parameters are also specific to British National Grid.

    fprintf('Northing = %.6f, Easting = %.6f\n', northing, easting);

    % Ellipsoid parameters (from Java code - appears to be OSGB36)
    a = 6377563.396;           % Semi-major axis
    fm1 = 299.3249646;         % Inverse flattening (1/f)
    f = 1.0 / fm1;             % Flattening
    e2 = 2*f - f*f;            % First eccentricity squared
    e = sqrt(e2);              % First eccentricity

    % Projection-specific parameters (British National Grid)
    projectionOriginLongitude = deg2rad(-2.0);
    projectionOriginLatitude = deg2rad(49.0);

    FE = 400000.0;             % False easting
    FN = -100000.0;            % False northing
    k0 = 0.9996012717;         % Scale factor at natural origin

    % Calculate n and powers of n
    n = f / (2 - f);
    dump_value(n, 'n');
    n2 = n^2;
    n3 = n^3;
    n4 = n^4;

    % Calculate B
    nB = a / (1 + n);
    dB = 1 + n2/4 + n4/64;
    B = nB * dB;
    dump_value(B, 'B');

    % Calculate h coefficients (non-dashed versions for M0 calculation)
    h1 = n/2 - (2/3)*n2 + (5/16)*n3 + (41/180)*n4;
    h2 = (13/48)*n2 - (3/5)*n3 + (557/1440)*n4;
    h3 = (61/240)*n3 - (103/140)*n4;
    h4 = (49561/161280)*n4;

    % Calculate h' (dashed) coefficients
    h1d = n/2 - (2/3)*n2 + (37/96)*n3 - (1/360)*n4;
    h2d = (1/48)*n2 - (1/15)*n3 + (437/1440)*n4;
    dump_value(h2d, 'h2 dash');
    h3d = (17/480)*n3 - (37/840)*n4;
    dump_value(h3d, 'h3 dash');
    h4d = (4397/161280)*n4;
    dump_value(h4d, 'h4 dash');

    % Calculate M0
    M0 = getM0(B, h1, h2, h3, h4, projectionOriginLatitude, e);

    % Calculate normalized coordinates
    ny = (easting - FE) / (B * k0);
    nE = (northing - FN + k0 * M0) / (B * k0);

    dump_value(k0, 'k0');
    dump_value(nE, 'nE');

    % Calculate e1d to e4d terms
    e1d = h1d * sin(2*nE) * cosh(2*ny);
    e2d = h2d * sin(4*nE) * cosh(4*ny);
    e3d = h3d * sin(6*nE) * cosh(6*ny);
    e4d = h4d * sin(8*nE) * cosh(8*ny);

    dump_value(e1d, 'e1d');
    dump_value(e2d, 'e2d');
    dump_value(e3d, 'e3d');
    dump_value(e4d, 'e4d');

    % Calculate n1d to n4d terms
    n1d = h1d * cos(2*nE) * sinh(2*ny);
    n2d = h2d * cos(4*nE) * sinh(4*ny);
    n3d = h3d * cos(6*nE) * sinh(6*ny);
    n4d = h4d * cos(8*nE) * sinh(8*ny);

    dump_value(n1d, 'n1d');
    dump_value(n2d, 'n2d');
    dump_value(n3d, 'n3d');
    dump_value(n4d, 'n4d');

    % Calculate n0d and e0d
    n0d = ny - n1d - n2d - n3d - n4d;
    e0d = nE - e1d - e2d - e3d - e4d;

    dump_value(n0d, 'n0d');
    dump_value(e0d, 'e0d');

    % Calculate conformal latitude Bd
    Bd = asin(sin(e0d) / cosh(n0d));
    dump_value(Bd, 'Bd');

    % Calculate initial Qd
    Qd = asinh(tan(Bd));

    % Iterative refinement (6 iterations as in Java code)
    IO = Qd;
    for i = 0:5
        diff = iter_function(e, IO);
        IO = Qd + diff;
        dump_value(IO, sprintf('QD%d', i));
    end

    % Also calculate using successive substitution (as in Java)
    Qdd = Qd + iter_function(e, Qd);
    Qddd = Qd + iter_function(e, Qdd);
    Qdddd = Qd + iter_function(e, Qddd);

    dump_value(Qdd, 'Qdd(M)   ');
    dump_value(Qddd, 'Qddd(M)  ');
    dump_value(Qdddd, 'Qdddd(M) ');

    % Calculate final latitude and longitude
    phi = atan(sinh(IO));
    lbda = asin(tanh(n0d) / cos(Bd));
    lbda = projectionOriginLongitude + lbda;

    % Convert to degrees
    lat = rad2deg(phi);
    lon = rad2deg(lbda);

    dump_value(lat, 'phi');
    dump_value(lon, 'lambda');
end

function result = iter_function(e, Qd)
    % Iteration function (translates iter() method from Java)
    eQd1 = tanh(Qd);
    eQd2 = e * eQd1;
    eQd3 = e * atanh(eQd2);
    result = eQd3;
end

function M0 = getM0(B, h1, h2, h3, h4, phi0, e)
    % Calculate M0 (translates getM0() method from Java)
    dump_value(phi0, 'phi0');

    if phi0 == 0
        M0 = 0;
    else
        % Calculate Q0
        dump_value(phi0, '    Psi0');
        Q0 = asinh(tan(phi0));

        Q01 = atanh(e * sin(phi0));
        Q01 = e * Q01;

        Q0 = Q0 - Q01;
        dump_value(e * e, '    e^2');
        dump_value(Q0, '    Q0');

        % Calculate B0
        B0 = sinh(Q0);
        B0 = atan(B0);
        dump_value(B0, '    b0');

        % Calculate Eo0
        Eo0 = asin(sin(B0));
        dump_value(Eo0, '    Eo0');

        % Calculate Eo terms
        Eo1 = h1 * sin(2 * Eo0);
        Eo2 = h2 * sin(4 * Eo0);
        Eo3 = h3 * sin(6 * Eo0);
        Eo4 = h4 * sin(8 * Eo0);

        dump_value(Eo1, '    Eo1');
        dump_value(Eo2, '    Eo2');
        dump_value(Eo3, '    Eo3');
        dump_value(Eo4, '    Eo4');

        % Calculate E0
        E0 = Eo0 + Eo1 + Eo2 + Eo3 + Eo4;

        dump_value(E0, '    E0');
        dump_value(B, '    B');

        % Calculate M0
        M0 = B * E0;
        dump_value(M0, '    M0');
    end
end

function dump_value(value, name)
    % Utility function to print values (translates dump() method from Java)
    fprintf('%s = %.15f\n', name, value);
end
