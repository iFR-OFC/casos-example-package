% Helper function to allow both double and CasADi symbolic variables
function mustBeNumericOrCasadi(x)
    if ~isa(x, 'double') && ~isa(x, 'casadi.SX')
        error('Input must be either double or casadi.SX.');
    end
end