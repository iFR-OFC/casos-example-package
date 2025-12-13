classdef PolySatellite < saero.Satellite
    %POLYSATELLITE Polynomial Approximation of Symbolic Satellite
    
    properties
        vi_B
        cosd = sym('cosd', 'real')
        significant_digits = 6;
        H1_fit_degrees;
        H2_fit_degrees;
    end
    
    methods
        function obj = PolySatellite(opts)
            %POLYSATELLITE
            arguments
                opts.vi_B ...
                        {saero.input_validation.mustBeNumericOrMatlabSym} ...
                        = sym('vi_B', [3,1], 'real')
                opts.inertia_B_B ...
                        {saero.input_validation.mustBeNumericOrMatlabSym} ...
                        = diag([1;2;2])
                opts.single_sided_panels = saero.Panels()
                opts.two_sided_panels = saero.Panels()
                opts.environment = saero.Environment()
                opts.H1_fit_degrees = [1,3];
                opts.H2_fit_degrees = [0, 2, 4];
            end
            obj.vi_B = opts.vi_B;
            obj.inertia_B_B = opts.inertia_B_B;
            obj.single_sided_panels = opts.single_sided_panels;
            obj.two_sided_panels = opts.two_sided_panels;
            obj.environment = opts.environment;
            obj.H1_fit_degrees = opts.H1_fit_degrees;
            obj.H2_fit_degrees = opts.H2_fit_degrees;

            % Validate Panels
            assert(obj.single_sided_panels.isSymOrNum(), ...
                "Single sided panels must be sym or double")
            assert(obj.two_sided_panels.isSymOrNum(), ...
                "Two sided panels must be sym or double")
        end

        function force_B = getAerodynamicForce(obj)
            warning(['Polynomial Model/fit currently only' ...
                ' considers two-sided panels'])
            % Polynomial fitting currently only uses 2 sided panels
            force_B = obj.getTwoSidedAerodynamicForce();
        end

        function torque_B = getAerodynamicTorque(obj)
            % Calculate Aerodynamic torque for each individual two-sided
            % panel
            force_B = obj.getAerodynamicForce();
            R_B = [obj.two_sided_panels.cop_pos]; % only use two-sided
            
            % Calculate cross product
            torque_B = cross(R_B,force_B);
        end

        function fig = plotFit(obj)
            H1 = matlabFunction(obj.getH1TwoSided);
            H2 = matlabFunction(obj.getH2TwoSided);
            H1Fit = matlabFunction(obj.getPolyFitH1TwoSided("fit_degrees", obj.H1_fit_degrees));
            H2Fit = matlabFunction(obj.getPolyFitH2TwoSided("fit_degrees", obj.H2_fit_degrees));

            fig = hfigure('Polynomial Fit vs original function'); clf; hold on;
            cd = linspace(-1,1,300);
            mlt.plot(cd, H1(cd));
            mlt.plot(cd, H2(cd));
            mlt.plot(cd, H1Fit(cd));
            mlt.plot(cd, H2Fit(cd));
            xlabel('$\cos \delta$')
            legend('$H1$','$H_{1}^P$','$H2$','$H_{2}^P$')
            hold off;
        end

        function fig = plotTorqueOverAlphaBeta(obj, opts)
            arguments
                obj 
                opts.figId  (1,:) char = 'Plot torque over alpha/beta'
                opts.alphaRange (1,:) double = linspace(-pi/2,pi/2,50);
                opts.betaRange (1,:) double = linspace(-pi/2,pi/2,50);
            end
            alpha = sym('a', 'real');
            beta = sym('b', 'real');

            obj.vi_B = ssmu.dcm.aeroToBody(alpha,beta)*[-1;0;0];

            torque = matlabFunction(obj.getTotalAerodynamicTorque);
            assert(nargin(torque)==2, ...
                'For plotting only vi_B can be symbolic')
            [Alpha, Beta] = meshgrid(opts.alphaRange, opts.betaRange);
            
            
            al = Alpha(:)';
            be = Beta(:)';
            tau = nan(3,size(al,2));
            for i=1:numel(Alpha)
                tau(:,i) = real(torque(al(i),be(i)));
            end
            
            
            % Plot
            fig = hfigure(opts.figId); hold on;
            tiledlayout(3,1);
            labels = {'\tau_x', '\tau_y', '\tau_z'};
            for i = 1:3
                ax = nexttile;
                surf(ax, Alpha, Beta, reshape(tau(i,:),size(Alpha)));
                xlabel('\alpha (rad)');
                ylabel('\beta (rad)');
                zlabel(labels{i});
                title(labels{i});
                % shading interp;
                colorbar;
            end
            hold off;
        end
        
        force_B = getTwoSidedAerodynamicForce(obj)
    end


    methods
        function H1Fit = getPolyFitH1TwoSided(obj, opts)
            % Fit analytic H1 function with matlab
            arguments
                obj 
                opts.fit_degrees (1,:) {mustBeInteger} = [1,3]
                opts.n_points (1,1) {mustBeInteger} = 500
            end
            % Get nonlinear two sided panel
            H1 = matlabFunction(obj.getH1TwoSided);

            assert(nargin(H1)==1,...
                "Number of inputs in H1 function mismatch")
            
            % Fit H1Poly
            H1Fit = obj.scalar_lsq_fit(H1, opts.fit_degrees, ...
                obj.cosd, [-1,1]);
        end

        function H2Fit = getPolyFitH2TwoSided(obj, opts)
            % Fit analytic H1 function with matlab
            arguments
                obj 
                opts.fit_degrees (1,:) {mustBeInteger} = [0,2,4]
                opts.n_points (1,1) {mustBeInteger} = 500
            end
            % Get nonlinear two sided panel
            H2 = matlabFunction(obj.getH2TwoSided);

            assert(nargin(H2)==1,...
                "Number of inputs in H1 function mismatch")
            
            % Fit H1Poly
            H2Fit = obj.scalar_lsq_fit(H2, opts.fit_degrees, ...
                obj.cosd, [-1,1]);
        end
    end
    methods (Static)
        function [poly, coeffs] = scalar_lsq_fit(fhandle, fit_degrees, ...
                indeterminate, range, opts)
            % Perform scalar lsq fit using symbolic toolbox
            arguments
                fhandle function_handle
                fit_degrees (1,:) {mustBeInteger}
                indeterminate (1,1) sym
                range (2,1) double
                opts.npoints = 500;
            end
            x = linspace(range(1), range(2), opts.npoints);
            z = real(fhandle(x)); % function values

            % Perform least-squares polynomial fitting with specified degrees
            H = zeros(length(x), length(fit_degrees));
            for i = 1:length(fit_degrees)
                H(:, i) = x.^fit_degrees(i);
            end

            % Solve least-squares problem
            coeffs = H \ z';
            
            % Get Polynomial as symexpr
            poly = 0;
            for i = 1:length(fit_degrees)
                poly = poly + coeffs(i)*indeterminate^(fit_degrees(i));
            end
        end
    end
end

