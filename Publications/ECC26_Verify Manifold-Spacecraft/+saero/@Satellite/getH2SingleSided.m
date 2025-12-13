function H2 = getH2SingleSided(obj)
    cosd = obj.cosd;
    si = obj.environment.si;
    
    % helper functions as function handles
    Gamma2 = @(x) 1./(sqrt(pi)).*exp(-x.^2) + x.*(1+erf(x));
                
    % H2 as symbolic expression
    H2 = 1./(si).*Gamma2(si.*cosd);
end