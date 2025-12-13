classdef (Abstract) Satellite
    %SATELLITE class to define a satellite geometry.
    
    properties (Abstract)
        vi_B (3,1) % Incoming velocity direction vector
        cosd (1,1) % Symbolic helper variable
    end

    properties
        inertia_B_B (3,3) = diag([1;2;2])
        single_sided_panels (1,1) saero.Panels
        two_sided_panels (1,1) saero.Panels
        environment (1,1) saero.Environment
    end
    
    methods
        H1 = getH1SingleSided(obj);
        H2 = getH2SingleSided(obj);
        H1 = getH1TwoSided(obj);
        H2 = getH2TwoSided(obj);

        function force_B_tot = getTotalAerodynamicForce(obj)
            force_B_tot = sum(obj.getAerodynamicForce, 2);
        end
        
        function torque_B_tot = getTotalAerodynamicTorque(obj)
            torque_B_tot = sum(obj.getAerodynamicTorque, 2);
        end
    end

    methods (Abstract)
        force_B = getAerodynamicForce(obj)
        torque_B = getAerodynamicTorque(obj)
    end
end

