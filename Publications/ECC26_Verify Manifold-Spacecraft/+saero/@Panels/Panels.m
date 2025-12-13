classdef Panels
    %PANELS Collection of (potentially symbolic) panels for a satellite in
    %free molecular flow
    
    properties
        isTwoSided logical = false
        normals (3,:)
        cop_pos (3,:)
        areas (1,:)
        wall_temperature (1,1)
    end
    
    methods
        function obj = Panels(normals, cop_pos, areas, opts)
            %PANELS Construct an instance of class Panels.
            arguments
                normals (3,:) {saero.input_validation.mustBeNumOrSym} = [];
                cop_pos (3,:) {saero.input_validation.mustBeNumOrSym, ...
                saero.input_validation.mustBeSameSize(normals, cop_pos)} = [];
                areas (1,:) {saero.input_validation.mustBeNumOrSym, ...
                saero.input_validation.mustBeSameLength(areas, normals, 2)} = [];
                opts.isTwoSided logical = false;
                opts.wall_temperature (1,1) ...
                    {saero.input_validation.mustBeNumOrSym} = 300;
            end
            obj.normals = normals;
            obj.cop_pos = cop_pos;
            obj.areas = areas;
            obj.isTwoSided = opts.isTwoSided;
            obj.wall_temperature = opts.wall_temperature;
        end
        
        function obj = addPanels(obj,newPanels)
            %ADDPANELS Adds set of panels to this set of panels
            arguments
                obj 
                newPanels (1,1) {mustBeA(newPanels, 'saero.Panels')}
            end

            % Validate that both panel sets have the same 'isTwoSided' 
            % property
            assert(obj.isTwoSided == newPanels.isTwoSided, ...
                "Cannot add panels: Both panel sets" ...
                + "must be either two-sided or single-sided.");
        
            % Check datatypes
            assert(isa(obj.normals, class(newPanels.normals)), ...
                'Panel normals must be of same type')
            assert(isa(obj.cop_pos, class(newPanels.cop_pos)), ...
                'Panel centres of pressure must be of same type')
            assert(isa(obj.areas, class(newPanels.areas)), ...
                'Panel centres of pressure must be of same type')

            obj.normals = [obj.normals, newPanels.normals];
            obj.cop_pos = [obj.cop_pos, newPanels.cop_pos];
            obj.areas = [obj.areas, newPanels.areas];
        end

        function isValid = isCasosOrNum(obj)
            %isCasosOrNum Checks if all properties are of type double or 
            % casos.PS
            allowedTypes = {'double', 'casos.PS'};
            
            % Collect property types
            propTypes = {class(obj.normals), class(obj.cop_pos), ...
                     class(obj.areas), class(obj.wall_temperature)};
            
            % Validate all types
            isValid = all(ismember(propTypes, allowedTypes));
        end

        function isValid = isSymOrNum(obj)
            %isCasosOrNum Checks if all properties are of type double or 
            % casos.PS
            allowedTypes = {'double', 'sym'};
            
            % Collect property types
            propTypes = {class(obj.normals), class(obj.cop_pos), ...
                     class(obj.areas), class(obj.wall_temperature)};
            
            % Validate all types
            isValid = all(ismember(propTypes, allowedTypes));
        end

        function isValid = isCasadiOrNum(obj)
            %isCasosOrNum Checks if all properties are of type double or 
            % casos.PS
            allowedTypes = {'double', 'casadi.SX'};
            
            % Collect property types
            propTypes = {class(obj.normals), class(obj.cop_pos), ...
                     class(obj.areas), class(obj.wall_temperature)};
            
            % Validate all types
            isValid = all(ismember(propTypes, allowedTypes));
        end
    end
end

