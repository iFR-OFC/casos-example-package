% Helper function to allow both double and CaSOS symbolic polynomials
function mustBeNumericOrCasosPS(x)
    if ~isa(x, 'double') && ~isa(x, 'casos.PS')
        error('Input must be either double or casos.PS');
    end
end