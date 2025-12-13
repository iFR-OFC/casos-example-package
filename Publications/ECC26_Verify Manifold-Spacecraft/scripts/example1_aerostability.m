%% Example 1: Aerostability of a cubesat with feathered geometry

%% User settings
plot_flag = true;

%% Kinematics and Dynamics
% State variables original coordinates
xbar        = casos.PD('xb', 6); % original coordinates
xbar_star1   = [1;zeros(5,1)];   % stable equilibrium in original coordinates
xbar_star2   = [-1;zeros(5,1)];  % unstable equilibrium in original coordinates

sbar       = xbar(1:3); % 2-axis attitude unit vector (negative wind direction)
ombar       = xbar(4:6); % rotational rates expressed in body coordinates

% Polynomial SOAR satellite
sat = saero.satellites.PolynomialSOAR();

% set reasonable environmental parameters in VLEO
sat.environment.updateParamsByAltitude("altitude_km", 300);

% shorthand: inertia matrix
J = sat.inertia_B_B;

% Aerodynamic torques using Sentman's method
torque.aeroFun = matlabFunction(sat.getTotalAerodynamicTorque());
torque.aero = torque.aeroFun(-sbar(2), -sbar(3)); % insert wind direction values

% Damping control law
Kd = 0.0160*J;
torque.control = -Kd*ombar;

% 2-axis attitude kinematics and dynamics
[sbar_dot, ombar_dot] = ssmu.models.dynamics.rigid_body_BI_2_axis( ...
    sbar, ombar, J, torque.aero+torque.control);

% System dynamics as symbolic expression
fbar = [sbar_dot; ombar_dot];

% manifold constraint in original coordinates
hbar = sbar'*sbar-1;

%% Shift and scale system dynamics
% x = S(xbar-xbar_star1)
S       = diag([ones(3,1);1/0.05.*ones(3,1)]);
x       = casos.PD('x', 6); % new state

% transform unstable equilibrium
xstar2 = mlt.transform.point.xbar_to_x(xbar_star2, xbar_star1, S);

% transform expressions
h = mlt.transform.expr.xbar_to_x(hbar, xbar, x, xbar_star1, S);

% transform dynamics (xdot = S*fbar(inv(S)*x+xbar_star1))
f = mlt.transform.system.xbar_to_x(fbar, xbar, x, xbar_star1, S);

%% Define SOS problem

% unknown polynomial decision variables
p = casos.PS.sym('p', monomials(x, 0:4)); % equality constraint 
V = casos.PS.sym('v', monomials(x, 1:1), 'gram'); % Lyap. candidate
sos.x = [p;V];
opts.Kx.lin = numel(p); % p is decision variable
opts.Kx.sos = numel(V); % V is SOS

% Enforce positive definite V
eps1 = 1e-5;
l1 = eps1*(x'*x);

% supply rate defined via equilibria
eps2 = 1e-5;
l2 = eps2*(x'*x)*((x-xstar2)'*(x-xstar2));

% sos constraints (SOS relaxation of Theorem 1)
sos.g = [
    V-l1;
    -nabla(V,x)*f+p'*h-l2
    ];
opts.Kc.sos = size(sos.g,1);

% Minimize Lyapunov coefficients (heuristic numerical improvement)
sos.f = dot(V,V);

% Turn off newton simplification
opts.newton_solver = [];

% Define Problem
Sos = casos.sossol('S', 'mosek', sos, opts);

% Solve
sol = Sos();

% Retrieve Lyapunov function solution and remove small (numerically 0)
% coefficients
Vsol = remove_coeffs(sol.x(2), 1e-10*max(sol.x(2).poly2basis));

%% Create Plot
if ~plot_flag
    return
end

% Lyapunov function in original coordinates
Vbar_sol = mlt.transform.expr.x_to_xbar(Vsol, x, xbar, xbar_star1, S);
Vbar_solFun = mlt.utils.toVecInFun(Vbar_sol.to_function);

% Dynamics as CasADi function
fbarFun = mlt.utils.toVecInFun(fbar.to_function);

% Random initial condition
% x0 = -1 + 2.*rand(6,1);
% x0(1:3) = x0(1:3)/norm(x0(1:3));
% x0(4:6) = 0.05.*x0(4:6); % scale rates to reasonable value

% Initial Condition used in paper:
x0 =[0.9969
   -0.0792
   -0.0006
   -0.0359
    0.0304
   -0.0105];

[t,y] = mlt.dynamics.ode45(fbarFun, x0, ...
    "duration__s", 800, "n_steps", 500);

% Evaluate Lyapunov function
Vbar_val = full(Vbar_solFun(y));

hfigure('Plot States and Lyapunov Function values (original coordinates');
clf;
hold on;
mlt.plot(t, y(1:3,:), LineStyle='--', LineWidth=.8);
mlt.plot(t, rad2deg(y(4:6,:)), LineWidth=.8);
mlt.plot(t, Vbar_val, 'k',LineWidth=1)
legend('$\bar{x}_1$','$\bar{x}_2$','$\bar{x}_3$', ...
    '$\bar{x}_4  \,[^{\circ}/s]$','$\bar{x}_5  \,[^{\circ}/s]$', ...
    '$\bar{x}_6  \,[^{\circ}/s]$', '$\bar{V}(\bar{x})$')
legend('Location', 'SouthEast')
hold off;



%% Analysis of unstable equilibria in transformed coordinates

% Linearize around unstable equilibrium
A_us = full(subs(nabla(f,x), x, xstar2));

% project Jacobian using into space tangential to H
ker_nabla_h_us = null(full(subs(nabla(h, x), x, xstar2)));
A_us_tang = ker_nabla_h_us'*A_us*ker_nabla_h_us;

% Eigenvalues of linerized system tangential to H
eig_us_tang = eig(A_us_tang)