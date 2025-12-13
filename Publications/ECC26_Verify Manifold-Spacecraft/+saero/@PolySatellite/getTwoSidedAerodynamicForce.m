function [force_B] = getTwoSidedAerodynamicForce(obj)
%SINGLESIDEDFORCE_B Calculate force of n single sided plates

% Get helper terms
H1 = obj.getPolyFitH1TwoSided("fit_degrees", obj.H1_fit_degrees);
H2 = obj.getPolyFitH2TwoSided("fit_degrees", obj.H2_fit_degrees);

% panels
panels = obj.two_sided_panels;
n = size(panels.normals,2);

% Calculate force for each panel
env = obj.environment;
N = panels.normals;
vi_B = repmat(obj.vi_B, 1, n);
A = panels.areas;

% Substitute H1 and H2 for each panel
H1vec = sym(zeros(size(A)));
H2vec = H1vec;
for i=1:n
    H1vec(i) = ssmu.subs(H1, obj.cosd, -N(:,i).'*vi_B(:,1));
    H2vec(i) = ssmu.subs(H2, obj.cosd, -N(:,i).'*vi_B(:,1));
end

% limited significant digits
H1vec = vpa(repmat(H1vec,3,1),obj.significant_digits); 
H2vec = vpa(repmat(H2vec,3,1),obj.significant_digits); 
 

p_B = env.density/2.*env.Vi.^2.*(H1vec.*N + H2vec.*vi_B);
force_B = repmat(A,3,1).*p_B; % Repmat required for Casadi Implementation
end