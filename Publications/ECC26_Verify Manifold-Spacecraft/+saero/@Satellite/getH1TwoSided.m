function H1 = getH1TwoSided(obj)
    H1single = obj.getH1SingleSided();
    cosd = obj.cosd;

    H1 = H1single - ssmu.subs(H1single, cosd, -cosd);
end