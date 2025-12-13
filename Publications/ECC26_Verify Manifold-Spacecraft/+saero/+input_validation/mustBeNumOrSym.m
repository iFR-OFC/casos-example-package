% Helper function to allows for numerical values, casadi.SX expressions or
% matlab symexpr
function mustBeNumOrSym(x)
    if ~isa(x, 'double') && ~isa(x, 'casadi.SX') ....
            && ~isa(x, 'sym') && ~isa(x, 'casos.PS')
        error('Input must be one of: double|casadi.SX|sym|casos.PS.');
    end
end