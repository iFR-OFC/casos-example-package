function mustBeSameLength(A, B, dim)
    if size(A, dim) ~= size(B, dim)
        error('Length of inputs A and B along dimension %d do not match. A: %d, B: %d', ...
              dim, size(A, dim), size(B, dim));
    end
end