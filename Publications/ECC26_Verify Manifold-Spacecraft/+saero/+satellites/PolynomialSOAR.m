classdef PolynomialSOAR < saero.PolySatellite
    %POLYNOMIALSOAR Polynomial version of SOAR satellite
    
    properties
    end
    
    methods
        function obj = PolynomialSOAR(g)
            %POLYNOMIALSOAR 
            arguments
                g.body_length (1,1) double = 0.3405; % body length
                g.body_width (1,1) double = 0.1; % body width and height
                g.com_to_rear (1,1) double = 0.161; % distance com to rear
        
                g.wing_length (1,1) double = 0.57; % wing: length
                g.wing_width (1,1) double = 0.06; % wing: width
                g.wing_distance_from_center (1,1) double = 0.335; % wing: distance from center
            end
            %SYMBOLICSOAR

            % Body geometry
            l = g.body_length;
            d = g.body_width;
            dr = g.com_to_rear;

            % wing geometry
            lw = g.wing_length;
            ww = g.wing_width;
            dw = g.wing_distance_from_center;

            % Build body
            body = saero.Panels(isTwoSided=false);
            body.normals = [  1,-1, 0, 0, 0, 0;
                            0, 0, 1, 0,-1, 0;
                            0, 0, 0,-1, 0, 1];
            body.areas = [d^2, d^2, d*l*ones(1,4)];

            body.cop_pos = [  0.366-dr, -dr,  l/2-dr, l/2-dr, l/2-dr, l/2-dr;
            0,        0,    d/2,    0,      -d/2,   0;
            0,        0,    0,      -d/2,   0,      d/2];

            % Two Sided parts (Wings)
            % bottom - right - top - left
            
            wings = saero.Panels(isTwoSided=true);

            % Normal vectors
            wings.normals = [  0, 0, 0, 0;
                        -1, 0, 1, 0;
                         0, 1, 0,-1];
            
            % wing areas
            wings.areas = lw*ww.*ones(1,4);
            
            wings.cop_pos = [ -dr+ww/2, -dr+ww/2, -dr+ww/2, -dr+ww/2;
                        0,        dw,        0,       -dw;
                        dw,        0,        -dw,     0];

            
            % Add panels to satellite
            obj.single_sided_panels = body;
            obj.two_sided_panels = wings;

            % Overwrite inertia value
            obj.inertia_B_B = diag([0.0288,0.0392,0.0392]);
        end
    end
end

