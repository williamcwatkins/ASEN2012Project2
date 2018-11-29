function [ProjectEquations] = Project2Equations(t,Initial)
    %Project2Equations defines the differential equations for the flight of
    % the bottle rocket
    %% Define the global variables
    global g Cd rhoAirAmb VolBottle PAmb gamma rhoWater  AThroat ...
        ABottle R CD PAbs x0 z0 Ls VolAirInitial mAirInitial ...
        TAirInitial InitialHeading Phase1Endt Phase2Endt

    %% Initial conditions for this timestep (not the FIRST timestep)
    x = Initial(1); % Current x position for the current timestep
    z = Initial(2); % Current z position for the current timestep
    vx = Initial(3); % Current x-component of velocity
    vz = Initial(4); % Current z-component of velocity
    mRocket = Initial(5); % Current mass of the whole rocket
    mAir = Initial(6); % Current mass of air in the rocket
    VolAir = Initial(7); % Current volume of the air in the rocket (increases to volume of bottle)

    %% Non-differential equations
    
    vMag = sqrt(vx ^ 2 + vz ^ 2); % Compute the speed of the rocket
    %disp(vMag);
    %disp('.');
    
    Drag = (rhoAirAmb / 2) * (vMag ^ 2) * CD * ABottle; % Equation 2, [N], 
    % Drag in Newtons on the bottle. Use Drag(num) to compute. 
    %disp(Drag);

    %% Calculate Pressures
    PAirThrust = 0;
    
    if VolAir >= VolBottle
        % Equation 13, Pressure and Temperature at end of water thrust phase
        PWaterThrustEnd = PAbs * (VolAirInitial / VolBottle) ^ gamma;
        TWaterThrustEnd = TAirInitial * (VolAirInitial / VolBottle) ^ (gamma - 1);
        % Equation 14 (rearranged), Pressure during air thrust phase at any point in time
        PAirThrust = PWaterThrustEnd * (mAir / mAirInitial) ^ gamma;
        % Equations 15, corresponding density and temperature with the
        % above pressure
        rhoAirThrust = mAir / VolBottle;
        TAirThrust = PAirThrust / (rhoAirThrust * R);
    end

    %% ALL the calculations!

    if VolAir < VolBottle % Phase 1, Water Thrust
        
        % Non-differential equations specific to Water Thrust Phase
        
        % Equation 3, Calculate pressure in [Pa] at current timestep
        PWaterThrust = PAbs * (VolAirInitial / VolAir) ^ gamma; 
        % Equation 5, Calculate thrust of the rocket in [N] at current timestep
        Thrust = 2 * Cd * AThroat * (PWaterThrust - PAmb);
        
        % Differential equations specific to Water Thrust Phase
        
        % Equation 9, Change in volume
        dVolAirdt = Cd * AThroat * sqrt((2 / rhoWater) * (PAbs * ((VolAirInitial / VolAir) ^ gamma) - PAmb));
        % Equation 10, Change in mass of rocket
        dmRocketdt = - Cd * rhoWater * AThroat * sqrt((2 * (PWaterThrust - PAmb)) / rhoWater);
        % For ODE45 to run, we also have to declare Equation 24, Change in air
        % mass, as 0
        dmAirdt = 0;

    elseif PAirThrust > PAmb % Phase 2, Air Thrust
        
        % Equation 16, Critical pressure in the bottle
        PCritical = PAirThrust * (2 / (gamma + 1)) ^ (gamma / (gamma - 1));

        if PCritical > PAmb % Calculations for a choked flow at the throat

            MExit = 1;
            % From Equation 18
            PExit = PCritical;
            % Equations 18, Exit temperature and density of the air
            TExit = (2 / (gamma + 1)) * TAirThrust;
            rhoExit = PExit / (R * TExit);
            % Equation 17, Exit velocity of the air
            vExit = sqrt(gamma * R * TExit);
            
        else % Calculations for a non-choked flow
            % Equation 19, Exit pressure is equal to ambient pressure
            PExit = PAmb;
            % Equation 20, Exit mach number
            MExit = sqrt((2 / (gamma - 1)) * (((PAirThrust / PAmb) ^ ((gamma - 1) / gamma)) - 1));
            % Equations 21, Exit Temperature and density
            TExit = TAirThrust / (1 + ((gamma - 1) / 2) * MExit ^ 2);
            rhoExit = PExit / (R * TExit);
            % Equation 22, Exit velocity of the air
            vExit = MExit * sqrt(gamma * R * TExit);
        end
        
        % Equation 24, Change in mass of air in rocket
        dmAirdt = -1* (Cd * rhoExit * AThroat * vExit);
        % Equation 23, Thrust developed by the rocket at current timesteo
        Thrust = (-1 * dmAirdt * vExit) + (PExit - PAmb) * AThroat;
        % Equation 26, Change in the total rocket mass
        dmRocketdt = dmAirdt;
        % To make ODE45 work, have to declare change in volume is 0
        dVolAirdt = 0;
    else % Phase 3, Ballistic trajectory
        % No thrust is produced
        Thrust = 0;
        % Volume of air is not changing
        dVolAirdt = 0;
        % Mass of rocket is constant
        dmRocketdt = 0;
        % Mass of air is constant
        dmAirdt = 0;
    end
    %% Velocity components
    % Get x- and y-components of the velocity
    dxdt = vx;
    dzdt = vz;
        
    %% Heading components
    % Check if the rocket is off the rail yet
    if sqrt((x - x0) ^ 2 + (z - z0) ^ 2) < Ls
        Hx = InitialHeading(1); % Initial x Heading
        Hz = InitialHeading(2); % Initial z Heading
    else % If the rocket has left the test stand
        Hx = dxdt / vMag; % x Heading
        Hz = dzdt / vMag; % z Heading
    end
    
    %% You were supposed to bring balance to the force!
    % x-component (multiplied by x-component of heading and divided by mass
    % to get acceleration in the x-direction)
    dvxdt = (Thrust - Drag) * (Hx / mRocket); 
    % y-component (multiplied by y-component of heading and divided by mass
    % to get acceleration in the y-direction)   
    dvzdt = (Thrust - Drag) * (Hz / mRocket) - g;

    %disp(Thrust)
    %disp('.')
    %disp(Drag)
    %disp(Hz)
    %disp(mRocket)
    
    %% Set our differential equations for ODE45
    
    ProjectEquations(1) = dxdt;
    ProjectEquations(2) = dzdt;
    ProjectEquations(3) = dvxdt;
    ProjectEquations(4) = dvzdt;
    ProjectEquations(5) = dmRocketdt;
    ProjectEquations(6) = dmAirdt;
    ProjectEquations(7) = dVolAirdt;
    ProjectEquations = ProjectEquations';
    
end