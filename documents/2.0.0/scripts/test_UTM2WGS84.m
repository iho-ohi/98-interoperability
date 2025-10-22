%% Test script for UTM2WGS84 function
% This tests the translated Java reverse() method

clear;
clc;

fprintf('=======================================================\n');
fprintf('Testing UTM2WGS84 (Translated from Proj.java)\n');
fprintf('=======================================================\n\n');

%% Test Case: Example coordinates
% Note: The Java implementation uses British National Grid parameters,
% not standard UTM parameters

fprintf('Test Case: British National Grid Coordinates\n');
fprintf('---------------------------------------------\n');

% Example British National Grid coordinates for a location
easting = 377178.0;    % Example easting
northing = 145794.0;   % Example northing

fprintf('\nInput Projected Coordinates:\n');
fprintf('  Easting:  %.6f m\n', easting);
fprintf('  Northing: %.6f m\n\n', northing);

% Convert to geographic coordinates
[lat, lon] = UTM2WGS84(easting, northing);

fprintf('\n=======================================================\n');
fprintf('FINAL RESULTS:\n');
fprintf('  Latitude:  %.10f°\n', lat);
fprintf('  Longitude: %.10f°\n', lon);
fprintf('=======================================================\n');
