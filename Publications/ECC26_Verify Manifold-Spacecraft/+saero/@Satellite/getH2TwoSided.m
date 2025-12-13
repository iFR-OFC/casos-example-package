function H2 = getH2TwoSided(obj)
    H2single = obj.getH2SingleSided();
    cosd = obj.cosd;

    H2 = H2single + ssmu.subs(H2single, cosd, -cosd);
end