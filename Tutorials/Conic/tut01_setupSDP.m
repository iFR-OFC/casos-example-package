%% -----------------------------------------------------------------------
%
% Short description: This tutorial describes how to setup a Semidefinite 
%                    Programming (SDP) problem in CaSoS. We guide you 
%                    through the process of defining variables, constraints,
%                    and solving the problem using CaSoS.
%
%  Problem: Minimize trace(C1*X1) + trace(C2*X2) subject to:
%           - trace(A1*X1) + trace(A2*X2) = b
%           - X2(1,2) >= -k  (since k = -3, this becomes X2(1,2) >= 3)
%           - X1(1,1)*[0,1; 1,3] + X2(2,2)*[3,1; 1,0] - I(2) >= 0
%           - X1 >= 0, X2 >= 0 (PSD constraints)
%
%   License: see License file of repository
%
% -----------------------------------------------------------------------

%% Step 1) Define matrix decision variables

% PSD matrix variables of different dimensions
X1 = casadi.SX.sym('X1', 3, 3);
X2 = casadi.SX.sym('X2', 4, 4);

%% Step 2) Define problem data (matrices and scalars)

% Cost matrices
C1 = sparse([1,3],[1,3], [1,6], 3, 3, 2);
C2 = sparse([1,1,2,2,3], [1,2,1,2,3], [1,-3,-3,2,1], 4, 4, 5);

% Constraint matrices
A1 = sparse([1,1,3,3],[1,3,1,3], [1,1,1,2], 3, 3, 4);
A2 = sparse([1,2,2,4], [2,1,2,4], [1,1,-1,-3], 4, 4, 4);

% Problem scalars
b = 23;
k = -3;

%% Step 3) Setup the problem struct

% Step 3.1: Define objective function
sdp.f = trace(C1*X1) + trace(C2*X2);

% Step 3.2: Define constraints
% Equality constraint: trace(A1*X1) + trace(A2*X2) - b = 0
g1 = trace(A1*X1)+trace(A2*X2)-b; 

% Inequality constraint: k - X2(1,2) >= 0  (i.e., X2(1,2) <= k)
g2 = k-(X2(1,2)+X2(2,1))/2;

% Matrix inequality (PSD constraint)
g3 = X1(1,1)*[0, 1; 1, 3] + X2(2,2)*[3, 1; 1, 0] - 1*eye(2);

% Collect all constraints
sdp.g = [g1; g2; g3(:)];

% Step 3.3: Define decision variables (vectorized)
sdp.x = [X1(:); X2(:)];

%% Step 4) Define cone structure

opts = struct;

% Decision variables cones: X1 and X2 are PSD matrices
opts.Kx = struct('psd', [3, 4]); 

% Constraint cones: 
% - g1 is linear (equality)
% - g2 is linear inequality (non-negative)
% - g3 is PSD constraint (2x2 matrix)
opts.Kc = struct('lin', 2, 'psd', 2);

%% Step 5) Generate a CaSoS solver instance

sdp_solver = 'mosek';

S = casos.sdpsol('S',        ...    % name of solver
                 sdp_solver, ...    % SDP solver
                 sdp,        ...    % problem structure
                 opts);             % options to solver (includes the cones definition)


%% Step 6) Call the solver to solve the SDP

fprintf('\n--- Solving SDP problem ---\n');
solve_start = tic;
sol = S('lbg', [0; 0],'ubg', [0; +inf]);
solve_time = toc(solve_start);

%% Step 7) Check the solution status

% Check with unified return status
if strcmp(S.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS')
    fprintf('  Problem solved successfully!\n');
elseif strcmp(S.stats.UNIFIED_RETURN_STATUS, 'SOLVER_RET_SUCCESS_REDUCED')
    fprintf('  Problem solved but with reduced accuracy\n');
else
    fprintf('  Solver failed: %s\n', S.stats.UNIFIED_RETURN_STATUS);
end

%% Step 8) Extract and display the solution

fprintf('\n--- Solution ---\n');

% Display objective value
if isfield(sol, 'f')
    fprintf('Optimal objective value: %.6f\n', full(sol.f));
end

% Reconstruct matrix solutions
n1 = 3; n2 = 4;
X1_sol = reshape(full(sol.x(1:n1^2)), n1, n1);
X2_sol = reshape(full(sol.x(n1^2+1:end)), n2, n2);

fprintf('\nX1 solution (3x3 PSD matrix):\n');
disp(X1_sol);

fprintf('X2 solution (4x4 PSD matrix):\n');
disp(X2_sol);


fprintf('\nConstraint values:\n');
fprintf('  g1 (equality): %.6f (should be 0)\n', full(sol.g(1)));
fprintf('  g2 (inequality): %.6f (should be >= 0)\n', full(sol.g(2)));

% Reshape g3 back to 2x2 matrix
g3_sol = reshape(full(sol.g(3:end)), 2, 2);
fprintf('  g3 (PSD constraint):\n');
disp(g3_sol);
fprintf('  Minimum eigenvalue of g3: %.6f\n', min(eig(g3_sol)));

