%% Example Usage of UTM and UPS to WGS84 Conversion Functions
% This script demonstrates how to use the utm_to_wgs84 and ups_north_to_wgs84 functions

clear;
clc;

fprintf('=== Coordinate Conversion Examples ===\n\n');

%% Example 1: UTM to WGS84 (UTM Zone 45S)
fprintf('Example 1: UTM Zone 45S to WGS84\n');
fprintf('----------------------------------\n');

% Input coordinates (example values)
E_utm = 500000;      % Easting in meters (center of zone)
N_utm = 5000000;     % Northing in meters
zone_number = 45;    % UTM Zone 45
is_north = false;    % South zone

% Convert to WGS84
[lat1, lon1] = utm_to_wgs84(E_utm, N_utm, zone_number, is_north);

fprintf('Input UTM Coordinates:\n');
fprintf('  Zone: %dS\n', zone_number);
fprintf('  Easting (E):  %.3f m\n', E_utm);
fprintf('  Northing (N): %.3f m\n', N_utm);
fprintf('\nOutput WGS84 Coordinates:\n');
fprintf('  Latitude:  %.8f°\n', lat1);
fprintf('  Longitude: %.8f°\n', lon1);
fprintf('\n');

%% Example 2: UTM to WGS84 (UTM Zone 18N - New York area)
fprintf('Example 2: UTM Zone 18N to WGS84 (New York area)\n');
fprintf('------------------------------------------------\n');

% Input coordinates (approximate for New York City)
E_utm2 = 583960;     % Easting in meters
N_utm2 = 4507523;    % Northing in meters
zone_number2 = 18;   % UTM Zone 18
is_north2 = true;    % North zone

% Convert to WGS84
[lat2, lon2] = utm_to_wgs84(E_utm2, N_utm2, zone_number2, is_north2);

fprintf('Input UTM Coordinates:\n');
fprintf('  Zone: %dN\n', zone_number2);
fprintf('  Easting (E):  %.3f m\n', E_utm2);
fprintf('  Northing (N): %.3f m\n', N_utm2);
fprintf('\nOutput WGS84 Coordinates:\n');
fprintf('  Latitude:  %.8f°\n', lat2);
fprintf('  Longitude: %.8f°\n', lon2);
fprintf('\n');

%% Example 3: UPS North to WGS84
fprintf('Example 3: UPS North to WGS84 (Arctic region)\n');
fprintf('---------------------------------------------\n');

% Input coordinates (example Arctic location)
E_ups = 2500000;     % Easting in meters
N_ups = 1500000;     % Northing in meters

% Convert to WGS84
[lat3, lon3] = ups_north_to_wgs84(E_ups, N_ups);

fprintf('Input UPS North Coordinates:\n');
fprintf('  Easting (E):  %.3f m\n', E_ups);
fprintf('  Northing (N): %.3f m\n', N_ups);
fprintf('\nOutput WGS84 Coordinates:\n');
fprintf('  Latitude:  %.8f°\n', lat3);
fprintf('  Longitude: %.8f°\n', lon3);
fprintf('\n');

%% Example 4: UPS North at the pole center
fprintf('Example 4: UPS North at origin (North Pole region)\n');
fprintf('--------------------------------------------------\n');

% Input coordinates (center of UPS projection)
E_ups2 = 2000000;    % Easting in meters (center)
N_ups2 = 2000000;    % Northing in meters (center)

% Convert to WGS84
[lat4, lon4] = ups_north_to_wgs84(E_ups2, N_ups2);

fprintf('Input UPS North Coordinates:\n');
fprintf('  Easting (E):  %.3f m\n', E_ups2);
fprintf('  Northing (N): %.3f m\n', N_ups2);
fprintf('\nOutput WGS84 Coordinates:\n');
fprintf('  Latitude:  %.8f°\n', lat4);
fprintf('  Longitude: %.8f°\n', lon4);
fprintf('\n');

fprintf('=== Conversion Complete ===\n');
