function [] = PlotParameters(t,Results)
%ThrustProfile outputs a graph of the thrust with respect to time
    %% Define the global variables
    global Cd VolBottle PAmb gamma AThroat R PAbs VolAirInitial mAirInitial ...
        TAirInitial
    stop = find(Results(:,2) < 0) - 1;
    t = t(1:stop);
    xvec = Results(1:stop,1);
    zvec = Results(1:stop,2);
    vxvec = Results(1:stop,3);
    vzvec = Results(1:stop,4);
    mRocketvec = Results(1:stop,5);
    mAirvec = Results(1:stop,6);
    VolAirvec = Results(1:stop,7);
    Thrust = zeros(1,length(t));
    Phase1Endt = 0;
    Phase2Endt = 0;
    
    for i = 1:length(t)
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

            if Phase1Endt == 0
                Phase1Endt = t(i-1);
            end
            
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
            % Ending time of Phase 2
            if Phase2Endt == 0
                Phase2Endt = t(i-1);
            end
            % No thrust is produced
            Thrust(i) = 0;

        end
    end
    
    %% Max positions
    xmax = max(xvec);
    zmax = max(zvec);
    ThrustMax = max(Thrust);    
    
    %% Set our limits
    xposstop = ceil(xmax + 2);
    zposstop = ceil(zmax + 2);
    ThrustStop = ceil(ThrustMax + 10);
    
    %% Create our vertical lines to delineate phase transitions
    Phase1Thrust = [0 ThrustStop];
    Phase1linet = [Phase1Endt Phase1Endt];
    
    Phase2Thrust = [0 ThrustStop];
    Phase2linet = [Phase2Endt Phase2Endt];
    
    % Create our plots
    figure(1);
    
    % First subplot, x and z position
    subplot(2, 1, 1)
    plot(xvec,zvec);
    title('Rocket Trajectory');
    xlim([0 xposstop]);
    xlabel('Distance [m]');
    ylim([0 zposstop]);
    ylabel('Height [m]');

    % Next subplot, thrust over time
    subplot(2, 1, 2)
    plot(t,Thrust)
    title('Thrust over time');
    xlim([0,0.5]);
    xlabel('Time [s]');
    ylim([0 ThrustStop]);
    ylabel('Thrust [N]');
    hold on
    plot(Phase1linet, Phase1Thrust);
    hold on
    plot(Phase2linet, Phase2Thrust);
    legend('Thrust profile','Phase 1 Transition Time','Phase 2 Transition Time');
    
    saveas(gcf,'Verification Case Rocket Trajectory.png');
    
end