close all;
clear;
%% Epoch Date and Decimal Year
epoch.Year = 2024;
epoch.Month = 10;
epoch.Day = 12;
epoch.Hour = 12;
epoch.Minute = 24;
epoch.Second = 45;
epoch.DecimalYear = decyear(epoch.Year, epoch.Month, epoch.Day, epoch.Hour, epoch.Minute, epoch.Second);
%% Orbital Elements
Semi_major_axis = 6921000;     % (m)
Inclination = 55;              % (deg)
RAAN = 0;                      % (deg)
Argument_of_latitude = 0;      % (deg)
simulator.EarthGConst = 3.986004418e14; % Earth's standard gravitational parameter in m^3/s^2
%% Attitude Kinematics & Dynamics
% Initial Euler Angles (rad)
attitude.Phi0 = 12 * pi /180;
attitude.Theta0 = -15 * pi /180;
attitude.Psi0 = 18 * pi /180;
% Initial Omega in Body (rad/s)
attitude.P0 = 0.1;
attitude.Q0 = 0.15;
attitude.R0 = 0.2;
% Moment of Inertia
satellite.MOI = [0.0393 -0.0044 -0.0100; -0.0044 0.0397 0.0096; -0.0100 0.0096 0.0118];
satellite.InvMOI = inv(satellite.MOI);
%% Simulation
visSampleTime = 0.1;
%% Magnetoqures Model
satellite.MTQ.DPMoment = 0.5; % Nominal Dipole
%% Power
satellite.MTQ.Power = [1.5 1.5 1.5]; % W
satellite.Gyro.NominalVoltage = 3.3;
satellite.Gyro.NominalPower = 0.5;
satellite.MM.NominalPower = 0.4;
satellite.MM.NominalVoltage = 3.3;

%% Sensor Specifications

% GPS Specs
satellite.GPS.Frequency = 1;                  %Hz
satellite.GPS.Position.Bias = 5;              % m
satellite.GPS.Position.NoiseSTDDEV = 2.5;     %m
satellite.GPS.Velocity.Bias = 0.1;            %m/s

% Gyro Specs
satellite.Gyro.Frequency = 1;                 %1Hz
satellite.Gyro.ScaleFactor = 0.5;             % precent (%)
satellite.Gyro.Range = [-150, 150] ;          % deg/s
satellite.Gyro.Bias = 7;                      % deg/hr
satellite.Gyro.AngularRandomWalk = 0.16;      % deg/sqrt(hr)
satellite.Gyro.RateRandomWalk = 6;            % deg/sqrt(hr^3)
satellite.Gyro.BiasInstability = 0.5;         % deg/s

% Magnetometer (MM) Specs
satellite.MM.Frequency = 1;                   %Hz
satellite.MM.Accuracy = 1;                    % precent (%)
satellite.MM.Sensitivity = 1;                 % mV/nT
satellite.MM.Range = [-1, 1].*2e5;            % nT
satellite.MM.ScaleFactorDrift = 100;          % ppm/˚C
satellite.MM.ZeroFieldDrift = 5;              % nT/˚C
satellite.MM.ZeroFieldBias = 300;             % nT
satellite.MM.NoiseVariance = [200; 100; 150]; % nT^2

%% Detumbling Specs

% k_Bdot = [10000 10000 10000];
J_MOI = diag(satellite.MOI);
k_Bdot = ((4 * pi)/(97*60)) * (1 + sin(55*pi/180)) * J_MOI;
MTQ_Swithces = [1 1 1];

%% Nadir Pointing Specs

%Kp = [10 10 10];
Kp1 = 10;
Kp2 = 10;
Kp3 = 10;
%Kd = [10 10 10];
Kd1 = 0;
Kd2 = 0;
Kd3 = 0;

%% Attitude Determination
NoiseStd = [0.002; 0.01];  % Sun sensor (0.002 rad), Magnetometer (0.01 rad)
signQuat = 1;              % Positive scalar convention for quaternion

