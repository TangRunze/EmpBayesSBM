function Xhat = asge(A,d)
% asge calculates the adjacency spectral graph embedding Xhat.
    [V, D] = eigs(A,d,'LA');
    Xhat = V*sqrt(D);
end