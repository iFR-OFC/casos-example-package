% Helper function to allows for numerical values, casadi.SX expressions or
% matlab symexpr
function mustBeNumericOrMatlabSym(x)
    if ~isa(x, 'double') && ~isa(x, 'sym')
        error('Input must be one of: double|sym');
    end
end