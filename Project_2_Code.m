%% ASEN 2012 Computational Methods Project 2
%{
William Watkins
Ponder Stine
%}

%% Clean out the workspace
clc
clear all
close all

%% Declare global variables
global g Cd rhoAirAmb VolBottle PAmb gamma rhoWater dThroat AThroat dBottle ...
    ABottle R mBottle CD PAbs Ls tSpan VolAirInitial TAirInitial mAirInitial ...
    VolWaterInitial PGageInitial mRocketInitial x0 z0 InitialHeading ...

%% Set all of our constants
g = 9.81; % [m/s^2], acceleration due to gravity
Cd = 0.8; % Discharge coefficient
rhoAirAmb = 0.961; % [kg/m^3], ambient air density
VolBottle = 0.002; % [m^3], volume of empty bottle
PAmb = 12.1 * 6894.7572931783; % [Pa], atmospheric pressure
gamma = 1.4; % ratio of specific heat for air
rhoWater = 1000; % [kg/m^3], density of water
dThroat = 0.021; % [m], diameter of throat
dBottle = 0.105; % [m], diameter of bottle
R = 287; % [J/(kg*K)], gas constant of air
mBottle = 0.15; % [kg], mass of empty two-liter bottle w/ cones and fins
CD = 0.5; % Drag coefficient
AThroat = pi * (dThroat / 2) ^ 2 ; % [m^2], Area of the throat
ABottle = pi * (dBottle / 2) ^ 2 ; % [m^2], Cross sectional area of the bottle
Ls = 0.5; % [m], length of test stand

%% Initial conditions for verification case

PGageInitial = 50 * 6894.7572931783; % [Pa], initial gage pressure of air in bottle
PAbs = PGageInitial + PAmb; % [Pa], absolute pressure of air in the bottle
VolWaterInitial = 0.001; % [m^3], inital volume of water inside bottle
VolAirInitial = VolBottle - VolWaterInitial; % [m^3], initial volume of air in bottle
TAirInitial = 300; %[K], initial temperature of air
mAirInitial = (PAbs * (VolAirInitial)) / (R * TAirInitial);
v0 = 0; % [m/s], initial velocity of rocket
vx = 0;
vz = 0;
thetainitial = 45; % [degrees], initial angle of rocket
InitialHeading = [cosd(thetainitial), sind(thetainitial)];
x0 = 0.0; % [m], initial horizontal distance
z0 = 0.25; % [m], initial vertical height
mRocketInitial = mBottle + (rhoWater * VolWaterInitial) + mAirInitial; % [kg], total mass of the bottle rocket system


%% Initial conditions for hitting 75 meters
%{
PGageInitial = 67 * 6894.7572931783; % [Pa], initial gage pressure of air in bottle
PAbs = PGageInitial + PAmb; % [Pa], absolute pressure of air in the bottle
VolWaterInitial = 0.001; % [m^3], inital volume of water inside bottle
VolAirInitial = VolBottle - VolWaterInitial; % [m^3], initial volume of air in bottle
TAirInitial = 300; %[K], initial temperature of air
mAirInitial = (PAbs * (VolAirInitial)) / (R * TAirInitial);
v0 = 0; % [m/s], initial velocity of rocket
vx = 0;
vz = 0;
thetainitial = 45; % [degrees], initial angle of rocket
InitialHeading = [cosd(thetainitial), sind(thetainitial)];
x0 = 0.0; % [m], initial horizontal distance
z0 = 0.25; % [m], initial vertical height
mRocketInitial = mBottle + (rhoWater * VolWaterInitial) + mAirInitial; % [kg], total mass of the bottle rocket system
%}

%% Computing our differentials

% Define all of our initial conditions
Initials = [x0 z0 vx vz mRocketInitial mAirInitial VolAirInitial];
% Initialize end times of each phase so we can plot them later
Phase1Endt = 0; 
Phase2Endt = 0;

% Time span,from 0 to 5 at .001 seconds timestep
tSpan = 0:0.001:5; % [s], integration time


% Call our ode45 function
% tSpanPhase1 is time vector, VPhase1 is Volume vector
% tPhase1End is time where V is volume of bottle, VPhase1End is volume of bottle, IndexPhase1End is index
[t, Results] = ode45('Project2Equations', tSpan, Initials);

%% Plotting our results

% Pulling our results
PlotParameters(t,Results);
stop = find(Results(:,2) < 0);
x = Results(1:stop-1,1);
z = Results(1:stop-1,2);
%plot(x,z)

maxx = floor(max(x) * 10^2)/10^2;
maxz = floor(max(z) * 10^2)/10^2;

fprintf('Maximum horizontal distance: %4.2f meters.\n', maxx);

fprintf('Maximum height: %4.2f meters.\n', maxz);