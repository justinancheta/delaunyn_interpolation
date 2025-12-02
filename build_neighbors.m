function N = build_neighbors(S)
% BUILD_NEIGHBORS Construct neighbor list for general-dimensional simplices.
%
% Input:
%   S : nSimplices × (d+1) matrix. Each row lists vertex indices of a simplex.
%
% Output:
%   N : nSimplices × (d+1) matrix.  N(i,j) is the index of the simplex
%       adjacent to simplex i across the face opposite vertex S(i,j).
%       0 indicates a boundary face (no neighbor).

    [nSimplices, d1] = size(S);
    d = d1 - 1;
    
    % Each simplex contributes (d+1) faces.
    numFaces = nSimplices * (d+1);

    % Data structure for storing faces
    FaceVerts = zeros(numFaces, d);      % each face has d vertices
    FaceID    = zeros(numFaces, 2);      % [simplexIndex, oppositeVertexIndex]

    % Step 1: list all faces
    k = 1;
    for s = 1:nSimplices
        verts = S(s,:);
        for j = 1:(d+1)
            % face = all vertices except the j-th one
            face = verts([1:j-1, j+1:end]);
            FaceVerts(k,:) = sort(face);          % sort so identical faces match
            FaceID(k,:)    = [s, j];
            k = k + 1;
        end
    end

    % Sort faces lexicographically so identical faces are adjacent
    [FaceVertsSorted, idx] = sortrows(FaceVerts);
    FaceID = FaceID(idx,:);

    % Step 2: match adjacent faces
    N = zeros(nSimplices, d+1);

    p = 1;
    while p <= numFaces
        q = p + 1;

        % Check if the next face matches (same vertices)
        if q <= numFaces && isequal(FaceVertsSorted(p,:), FaceVertsSorted(q,:))
            s1 = FaceID(p,1);   j1 = FaceID(p,2);
            s2 = FaceID(q,1);   j2 = FaceID(q,2);

            % Assign neighbors
            N(s1,j1) = s2;
            N(s2,j2) = s1;

            p = p + 2; % skip matched pair
        else
            p = p + 1; % boundary face
        end
    end

end
