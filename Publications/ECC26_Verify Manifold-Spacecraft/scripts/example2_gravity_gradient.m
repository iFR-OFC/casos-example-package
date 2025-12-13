%% Example 2: Stabilization of a cubesat in circular orbit

%% User settings
plot_flag = true;

%% Kinematics and Dynamics
% State variables original coordinates
xbar        = casos.PD('xb', 7); % original coordinates
xbar_star1   = [1;zeros(6,1)];   % stable equilibrium in original coordinates
xbar_star2   = [-1;zeros(6,1)];  % unstable equilibrium in original coordinates

qbar       = xbar(1:4);  % Attitude quaternion B w.r.t. O
ombar       = xbar(5:7); % rotational rates expressed in body coordinates

% inertia matrix
J = diag([0.0288, 0.0392, 0.0392]);

% Direction Cosine Matrix
dcm = ssmu.quat.toDCM(qbar);

% Constants
constants.mu_earth = 3.986004418e14;  % m^3/s^2
constants.R_earth  = 6.378e6;         % m
constants.altitude_m = 500e3;         % m

% Orbital rate
omega0 = sqrt(constants.mu_earth / ...
    (constants.R_earth + constants.altitude_m)^3);

% Control law based on linearization from
% [Wie, 1989: Quaternion Feedback Regulator for Spacecraft Eigenaxis 
% Rotation - Section IV]
ctrl.t_s    = 50;                       % settling time
ctrl.Kd     = 4/ctrl.t_s*J;             % D-Gain
ctrl.Kp     = 2*2*(2/ctrl.t_s)^2*J;     % P-Gain
torque.control = -ctrl.Kd*ombar - ctrl.Kp*qbar(2:end);

% Gravity gradient torque
torque.gravity_gradient = ssmu.models.torque.gravity_gradient_circular( ...
    dcm(:,3), J, omega0);

[qbar_dot, ombar_dot] = ssmu.models.dynamics.rigid_body_BT_circular(qbar, ombar, J, ...
    torque.gravity_gradient+torque.control, omega0);

% Dynamics in original coordinates
fbar   = [qbar_dot; ombar_dot];

% manifold constraint in original coordinates
hbar = qbar'*qbar-1;


%% Shift and scale system dynamics
% x = S(xbar-xbar_star1)
S       = diag([ones(4,1);15.*ones(3,1)]);
x       = casos.PD('x', 7); % new state

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

% Enforce positive definiteness on V
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

% Turn off newton simplify
opts.newton_solver = [];


Sos = casos.sossol('S', 'mosek', sos, opts);

sol = Sos();

% Retrieve Lyapunov function solution and remove small (numerically 0)
% coefficients
Vsol = remove_coeffs(sol.x(2), 1e-3*max(sol.x(2).poly2basis));


%% Create Plot
if ~plot_flag
    return
end


% Lyapunov function in original coordinates
Vbar_sol = mlt.transform.expr.x_to_xbar(Vsol, x, xbar, xbar_star1, S);
Vbar_solFun = mlt.utils.toVecInFun(Vbar_sol.to_function);

% Dynamics as CasADi function
fbarFun = mlt.utils.toVecInFun(fbar.to_function);

% Random Initial condition
% x0 = -1 + 2.*rand(7,1);
% x0(1:4) = [-.93;0.1;0;0];
% x0(1:4) = x0(1:4)/norm(x0(1:4)); % Quaternion unit constraint
% x0(5:7) = 0.05.*x0(5:7); % scale rates to reasonable value

% Initial Condition from paper plot
x0 =[-0.9943
    0.1069
         0
         0
    0.0377
    0.0043
    0.0363];

[t,y] = mlt.dynamics.ode45(fbarFun, x0, ...
    "duration__s", 250, "n_steps", 700);

% Evaluate Lyapunov function
Vbar_val = full(Vbar_solFun(y));

hfigure('Plot States and Lyapunov Function values (original coordinates');
clf;
hold on;
mlt.plot(t, y(1:4,:), LineStyle='--', LineWidth=.8);
mlt.plot(t, rad2deg(y(5:7,:)), LineWidth=.8);
mlt.plot(t, Vbar_val, 'k',LineWidth=1)
legend('$\bar{x}_1$','$\bar{x}_2$','$\bar{x}_3$', ...
    '$\bar{x}_4  \,[^{\circ}/s]$','$\bar{x}_5$', ...
    '$\bar{x}_6  \,[^{\circ}/s]$','$\bar{x}_7  \,[^{\circ}/s]$', ...
    '$V(\bar{x})$', ...
    'Location', 'eastoutside')
hold off;

%% Analysis of unstable equilibria in transformed coordinates

% Linearize around unstable equilibrium
A_us = full(subs(nabla(f,x), x, xstar2));

% project Jacobian using into space tangential to H
ker_nabla_h_us = null(full(subs(nabla(h, x), x, xstar2)));
A_us_tang = ker_nabla_h_us'*A_us*ker_nabla_h_us;

% Eigenvalues of linerized system tangential to H
eig_us_tang = eig(A_us_tang)