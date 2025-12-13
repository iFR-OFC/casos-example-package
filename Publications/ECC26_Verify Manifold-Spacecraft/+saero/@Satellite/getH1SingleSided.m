function H1 = getH1SingleSided(obj)
    % Get environment params
    env = obj.environment;
    kB = env.kB;
    mT = env.mT;
    Vi = env.Vi;
    alphaE = env.alpha_E;
    si = env.si;
    
    cosd = obj.cosd;
    Tw = obj.single_sided_panels.wall_temperature;
    
    % Temperature ratio (1xn)
    enum = si.*cosd.*ssmu.erfc(-si.*cosd); % (1xn)
    denom = 1./sqrt(pi).*exp(-si.^2*cosd.^2) + ....
        si.*cosd.*ssmu.erfc(-si.*cosd);
    Trat = alphaE.*(2.*kB.*Tw)./(mT.*Vi.^2).*si^2 + ...
                    (1-alphaE).*( ...
                        1 + si.^2./2 + 1/4.* enum./denom...
                    );
    
    % helper functions as function handles
    Gamma1 = @(x) x./(sqrt(pi)).*exp(-x.^2) ...
        + (1/2 + x.^2).*(1+erf(x));
    Gamma2 = @(x) 1./(sqrt(pi)).*exp(-x.^2) + x.*(1+erf(x));
    
    % H1 as 1x1 symexpr
    H1 = 1./(si.^2).*(-Gamma1(si.*cosd) ...
        - sqrt(pi)/2.*sqrt(Trat).*Gamma2(si.*cosd) + ...
        si.*Gamma2(si.*cosd).*cosd);
end