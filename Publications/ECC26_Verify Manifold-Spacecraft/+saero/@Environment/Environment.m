classdef Environment < handle
    %ENVIRONMENT Class used for aerodynamic modelling
    
    properties
        density (1,1) {saero.input_validation.mustBeNumOrSym}
        Vi (1,1) {saero.input_validation.mustBeNumOrSym}
        alpha_E (1,1) {saero.input_validation.mustBeNumOrSym}
        si (1,1) {saero.input_validation.mustBeNumOrSym}
    end

    properties (Constant)
        mu_Earth = 3.986e14;
        r_Earth = 6.378e6;
        kB = 1.380649e-23;  % Boltzmann Constant
        mT = 16 * 1.6605390689252e-27;
    end
    
    methods
        function obj = Environment(opts)
            %ENVIRONMENT Construct an instance of class Environment
            arguments
                opts.density = 6e-11;
                opts.Vi = 7800;
                opts.alpha_E = 0.95;
                opts.si = 7;
            end
            obj.density = opts.density;
            obj.Vi = opts.Vi;
            obj.alpha_E = opts.alpha_E;
            obj.si = opts.si;

            obj = obj.updateParamsByAltitude("altitude_km", 200);
        end

        function obj = updateParamsByAltitude(obj, opts)
            arguments
                obj
                opts.altitude_km (1,1) double = 300;
                opts.lat_deg (1,1) double = 0;
                opts.lon_deg (1,1) double = 0;
                opts.year (1,1) {mustBeInteger} = 2024;
                opts.dayOfYear (1,1) {mustBeInteger} = 150;
                opts.UTseconds (1,1) {mustBeInteger} = 0;
            end
            [T, R] = atmosnrlmsise00( ...
                1000*opts.altitude_km, ...
                opts.lat_deg, opts.lon_deg, ...
                opts.year, opts.dayOfYear, opts.UTseconds);

            Ti = T(2);
            obj.density = R(6);

            cmi = sqrt(2*obj.kB*Ti/obj.mT);
            obj.Vi = sqrt(obj.mu_Earth / (obj.r_Earth + opts.altitude_km*1000));
            obj.si = obj.Vi / cmi;
        end

        function isValid = isCasosOrNum(obj)
            %isCasosOrNum Checks if all properties are of type double or 
            % casos.PS according to requirements
            allowedTypes = {'double', 'casos.PS'};
            
            % Collect property types
            propTypes = {class(obj.density), class(obj.Vi)};
            
            % Validate all types
            isValid = all(ismember(propTypes, allowedTypes));
        end
    end
end

