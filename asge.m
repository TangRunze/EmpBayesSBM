function Xhat = asge(A, dimLatentPosition)
% asge calculates the adjacency spectral graph embedding Xhat.
    [V, D] = eigs(A, dimLatentPosition, 'LA');
    Xhat = V*sqrt(D);
end