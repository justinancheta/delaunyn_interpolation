function lambda = barycentric_nd(q, Vsimplex)
% BARYCENTRIC_ND  Compute normalized barycentric coordinates of q
% inside a d-simplex defined by Vsimplex.
%
% Input:
%   q         : 1×d point
%   Vsimplex  : (d+1)×d simplex vertex matrix
%
% Output:
%   lambda    : 1×(d+1) barycentric coordinates (sum = 1)
%
% Compatible with old MATLAB versions (pre-R2016) using bsxfun.

    % Dimensions
    d = size(Vsimplex, 2);

    % -----------------------------------------------------------
    % Build matrix A = [v1 - vd+1, v2 - vd+1, ..., vd - vd+1]
    % -----------------------------------------------------------
    A = zeros(d, d);
    base = Vsimplex(d+1, :);     % last vertex v_{d+1}

    for i = 1:d
        A(:, i) = (Vsimplex(i, :) - base).';  % column vector
    end

    % RHS = q - v_{d+1}
    rhs = (q - base).';

    % Solve A * x = rhs   (x are the first d barycentric coordinates)
    % Robust solution with fallback
    if rcond(A) < 1e-12
        % Degenerate simplex ? return uniform weights
        lambda = ones(1, d+1) / (d+1);
        return;
    end

    x = A \ rhs;      % d×1

    % Last barycentric coordinate
    lambda = zeros(1, d+1);
    lambda(1:d) = x';
    lambda(d+1) = 1 - sum(x);

    % -----------------------------------------------------------
    % Normalize barycentric coordinates (important for numerics)
    % -----------------------------------------------------------
    s = sum(lambda);
    if abs(s) > 1e-14
        lambda = lambda / s;
    end

end
