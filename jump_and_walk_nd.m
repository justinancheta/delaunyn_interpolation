function simplexID = jump_and_walk_nd(V, S, N, q, startSimplex)
% V: nVertices × d   (coordinates)
% S: nSimplices × (d+1)
% N: nSimplices × (d+1), neighbors opposite each vertex
% q: 1 × d (query point)
% startSimplex: initial guess for the simplex (random or heuristic)
%
% Returns simplexID or 0 if outside.

    if nargin < 5
        startSimplex = randi(size(S,1));
    end
    simplexID = 0; % fail-safe
    current = startSimplex;
    MAX_STEPS = 2000;

    for step = 1:MAX_STEPS

        % Vertex coordinates of the current simplex
        verts = S(current,:);
        Vsimplex = V(verts,:);

        % Barycentric coordinates
        lambda = barycentric_nd(q, Vsimplex);

        % Check if inside
        if all(lambda >= -1e-12)
            simplexID = current;
            return;
        end

        % Find most negative barycentric coordinate
        [~, j] = min(lambda);

        next = N(current, j);

        if next == 0
            simplexID = 0; % point outside mesh
            return;
        end

        current = next;
    end

    
end
