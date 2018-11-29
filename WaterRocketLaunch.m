function[time_data,height_data,velocity_data] = WaterRocketLaunch(vol_bottle, percent_water, p, C_d, mode, bottle_dia)
% Water Bottle Rocket Project
% MAE 250, North Carolina State University
% Sean Murray, Conner Grey, James Einwaechter
% 12/12/2017
% Description:
%   Takes initial variables of a water bottle rocket and models it's height
%   and velocity over time
% Inputs: 
%   p: Rocket Pressure, (N/m^2)
%   percent_water: (volume of water in rocket)/(volume of rocket)  
%   vol_bottle: Volume of rocket, (m^3)
%   mode: 0 -model boosted phase
%         1 -model boosted and free flight phases
%         2 -model to peak altidude
%   bottle_dia: Diamater of Water Bottle Rocket (m)
% Outputs:
%   time_data: matrix with time of each loop iteration
%   height_data: matrix with height of each loop iteration
%   velocity_data: matrix with velocity of each loop iteration
% Initilize Variables, SI Units
t = 0;  
% Time, (s)
v = 0;  
% Velocity, (m/s)
h = 0;  
% Height, (m)
vol_water = vol_bottle*percent_water; 
% Volume of Water, (m^3)
vol_air = vol_bottle - vol_water; 
% Volume of Air, (m^3)
delta_t = 0; 
% Change in Time, (s)
delta_vol = 1e-5; 
% Change in Volume, (m^3)
time_data = [];
height_data = [];
velocity_data = [];
% Initilize Constants, SI Units
p_atm = 1.0018e5;   
% Pressure, (N/m^2)
g = 9.8064; 
% Gravity, (m/s^2)
rho_water = 1000;   
% Density of Water, (kg/m^3)
rho_air = 1.2137; 
% Density of Air, (kg/m^3)
rocket_mass = .074; 
% Mass of Rocket, (kg) 
d_hole = .0215; 
% Diameter of Hole, (m)
d_rocket = bottle_dia; 
% Diameter of Rocket, (m)
initial_vol_water = vol_water; 
% Initial Volume of Air (m^3))
% Loop Until Burnout
while(vol_water >= 0)
u = sqrt(2*p/rho_water)*sqrt(((vol_water/(initial_vol_water)^1.4) - p_atm/p)/(1-((d_hole^2)/(d_rocket^2)))^2);
delta_t = delta_vol/(u*pi/4*d_hole^2);
t = t + delta_t;
vol_air = vol_air + delta_vol;
vol_water = vol_water -delta_vol;
delta_v = -g*delta_t+(rho_water*delta_vol)/(rho_water*vol_water+rocket_mass)*u;
v = v + delta_v;
h = h + v*delta_t;
time_data = [time_data,real(t)];
height_data = [height_data,real(h)];
velocity_data = [velocity_data,real(v)];
end
if mode == 1
% Loop Until Landing
delta_t = .001;
while(h >= 0)
if v > 0    
% Rocket Climing
delta_v = (-g - .5*rho_air*v^2*(pi*(d_rocket^2)/4)*C_d/rocket_mass)*delta_t;
else
% Rocket Falling
delta_v = (-g+.5*rho_air*v^2*(pi*(d_rocket^2)/4)*C_d/rocket_mass)*delta_t;
end
v = v + delta_v;
h = h + v*delta_t;
t = t + delta_t;
time_data = [time_data,real(t)];
height_data = [height_data,real(h)];
velocity_data = [velocity_data,real(v)];
end
elseif mode == 2
% Loop Until Max Height
delta_t = .001;
while(v >= 0)
if v > 0    
% Rocket Climing
delta_v = (-g-.5*rho_air*v^2*(pi*(d_rocket^2)/4)*C_d/rocket_mass)*delta_t;
else
% Rocket Falling
delta_v = (-g+.5*rho_air*v^2*(pi*(d_rocket^2)/4)*C_d/rocket_mass)*delta_t;
end
v = v + delta_v;
h = h + v*delta_t;
t = t + delta_t;
time_data = [time_data,real(t)];
height_data = [height_data,real(h)];
velocity_data = [velocity_data,real(v)];
end
end
end