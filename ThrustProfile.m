function [] = PlotParameters(t,Results)
%ThrustProfile outputs a graph of the thrust with respect to time
    %% Define the global variables
    global g Cd rhoAirAmb VolBottle PAmb gamma rhoWater dThroat AThroat dBottle ...
        ABottle R mBottle CD PGageInitial PAbs v0 x0 z0 Ls VolAirInitial mAirInitial mRocketInitial ...
        TAirInitial InitialHeading PressureTotal ThrustTotal DragTotal i
    Thrust = zeros(1,length(t));
    for i = 1:length(Results)
        x = Results(i,1); % Current x position for the current timestep
        z = Results(i,2); % Current z position for the current timestep
        vx = Results(i,3); % Current x-component of velocity
        vz = Results(i,4); % Current z-component of velocity
        mRocket = Results(i,5); % Current mass of the whole rocket
        mAir = Results(i,6); % Current mass of air in the rocket
        VolAir = Results(i,7); % Current volume of the air in the rocket (increases to volume of bottle)
        
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
            Thrust(i) = 2 * Cd * AThroat * (PWaterThrust - PAmb);

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
            dmAirdt = - (Cd * rhoExit * AThroat * vExit);
            % Equation 23, Thrust developed by the rocket at current timesteo
            Thrust(i) = - (dmAirdt * vExit) + (PExit - PAmb) * AThroat;
            % Equation 26, Change in the total rocket mass

        else % Phase 3, Ballistic trajectory
            % No thrust is produced
            Thrust(i) = 0;

        end
    end
    %% Plot everything
    plot(t,Thrust)
    xlim([0,0.5])
end

